#!/usr/bin/env python
"""ANM mechanics features (protocol sections 10-12).

Frozen implementation:
- C-alpha coarse-graining, cutoff 10 A, spring constant 1.
- 3N x 3N Hessian; remove 6 rigid-body zero modes; full dense pseudo-inverse of the
  non-rigid subspace (eigendecomposition, eigenvalues > tol).
- R_ij = Frobenius norm of the 3x3 block (i,j) of H+.
- DFI_i = sum_{j != i} R_ij; %DFI = percentile rank within core residues.
- Dispersion: p_ij = R_ij / sum_{k != i} R_k; D_i = -sum p_ij log p_ij / log(N-1);
  percentile rank within core residues.
- core_residue_mask: trim from each terminus until first residue with >= 3 contacts (10 A).
"""
import numpy as np

CUTOFF = 10.0
TOL = 1e-8

def anm_hessian(ca, cutoff=CUTOFF):
    n = len(ca)
    d = np.linalg.norm(ca[:, None] - ca[None, :], axis=2)
    H = np.zeros((3 * n, 3 * n), dtype=np.float64)
    for i in range(n):
        js = np.where((d[i] <= cutoff) & (d[i] > 0))[0]
        for j in js:
            dij = d[i][j]
            diff = (ca[i] - ca[j]).astype(np.float64)
            block = np.outer(diff, diff) / dij ** 2
            H[3*i:3*i+3, 3*j:3*j+3] = -block
        diag = np.zeros((3, 3))
        for j in js:
            dij = d[i][j]
            diff = (ca[i] - ca[j]).astype(np.float64)
            diag += np.outer(diff, diff) / dij ** 2
        H[3*i:3*i+3, 3*i:3*i+3] = diag
    return H, d

def response_matrix(ca, cutoff=CUTOFF):
    """Return R (N,N) with R_ij = Frobenius norm of 3x3 block (i,j) of H+ (non-rigid)."""
    n = len(ca)
    H, d = anm_hessian(ca, cutoff)
    # symmetry QC
    sym_err = float(np.abs(H - H.T).max())
    lam, V = np.linalg.eigh(H)
    # rigid-body modes: eigenvalues ~ 0 (6 expected)
    zero = lam <= TOL
    n_zero = int(zero.sum())
    nz = ~zero
    W = V[:, nz] / np.sqrt(lam[nz])[None, :]  # (3N, m)
    # A_i = W_i W_i^T (3x3); R_ij = sqrt(trace(A_i A_j))
    A = W.reshape(n, 3, -1)
    # trace(A_i A_j) = sum_ab A_i[a,b] A_j[a,b]
    R2 = np.einsum("iab,jab->ij", A, A)
    R = np.sqrt(np.maximum(R2, 0))
    # QC: H H+ H ~ H on the non-rigid subspace
    Hplus = (V[:, nz] / lam[nz]) @ V[:, nz].T
    proj_err = float(np.abs(H @ Hplus @ H - H).max())
    qc = dict(sym_err=sym_err, n_zero_modes=n_zero, min_nonzero_eig=float(lam[nz].min()),
              max_nonzero_eig=float(lam[nz].max()), proj_err=proj_err, n_res=n)
    return R, qc

def core_mask(contacts, min_contacts=3):
    """Trim from each terminus until first residue with >= min_contacts contacts."""
    n = len(contacts)
    lo = 0
    while lo < n and contacts[lo] < min_contacts:
        lo += 1
    hi = n - 1
    while hi > lo and contacts[hi] < min_contacts:
        hi -= 1
    mask = np.zeros(n, dtype=bool)
    mask[lo:hi+1] = True
    return mask, lo, hi

def mechanics_features(ca, contacts, cutoff=CUTOFF):
    """Return dict with raw DFI, dispersion, core mask, QC."""
    R, qc = response_matrix(ca, cutoff)
    n = len(ca)
    Roff = R.copy()
    np.fill_diagonal(Roff, 0)
    dfi = Roff.sum(1)
    rowsum = Roff.sum(1)
    disp = np.zeros(n)
    valid = rowsum > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        p = Roff[valid] / rowsum[valid, None]
        logp = np.where(p > 0, p * np.log(p), 0.0)
        disp[valid] = -logp.sum(1) / np.log(n - 1)
    mask, lo, hi = core_mask(contacts)
    return dict(dfi=dfi, dispersion=disp, core_mask=mask, trim_lo=lo, trim_hi=hi, qc=qc, R=R)

def percentile_rank(x):
    """Percentile rank (0-100) with average-rank ties, as PercentileRank."""
    from scipy.stats import rankdata
    ranks = rankdata(x, method="average")
    return 100.0 * (ranks - 0.5) / len(x)
