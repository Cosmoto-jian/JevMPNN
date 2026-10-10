#!/usr/bin/env python
"""Static residue features (protocol section 8, feature_dictionary.tsv).

Implementation (frozen):
- SASA/RSA: Shrake-Rupley (numpy, 92 fibonacci points/atom, probe 1.4 A, element vdW radii),
  per-residue absolute SASA / theoretical maximum (Tien et al. 2013).
- Depth: accessible surface points of all heavy atoms; residue depth = mean over residue heavy
  atoms of distance to nearest accessible surface point.
- WCN: sum 1/d^2 over C-alpha within 10 A excluding self.
- SS: mkdssp 4.2.2, 3-state H/E/C.
- Radial position, terminal distance, 10 A neighborhood means.
"""
import os, subprocess, tempfile
import numpy as np
from scipy.spatial import cKDTree
import gemmi

MAX_ASA = {  # Tien et al. 2013
    "ALA":129,"ARG":274,"ASN":195,"ASP":193,"CYS":167,"GLN":225,"GLU":223,"GLY":104,"HIS":224,
    "ILE":197,"LEU":201,"LYS":236,"MET":224,"PHE":240,"PRO":159,"SER":155,"THR":172,"TRP":285,
    "TYR":263,"VAL":174,"MSE":224,"SEC":184,"PYL":241,"ASX":194,"GLX":224,"UNK":200,
}
AA20 = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO",
        "SER","THR","TRP","TYR","VAL"]
VDW = {"C":1.88,"N":1.64,"O":1.46,"S":1.77,"P":1.90,"SE":1.90}  # ProtOr radii
PROBE = 1.4
N_POINTS = 92

def read_structure(path):
    st = gemmi.read_structure(path)
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.remove_waters()
    return st

def chain_residues(st, chains):
    recs = []
    for chname in chains:
        ch = st[0][chname]
        for r in ch:
            info = gemmi.find_tabulated_residue(r.name)
            if info is not None and info.is_amino_acid():
                if r.find(" ", gemmi.IfSeqid.ANY) if False else ("CA" not in r):
                    continue  # residues without CA (rare partial records) are skipped
                recs.append((chname, r))
    return recs

def _fibonacci_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], 1)

def shrake_rupley(coords, elements, n_points=N_POINTS, probe=PROBE):
    """Per-atom SASA. coords (N,3) float32, elements list of element symbols."""
    radii = np.array([VDW.get(e.upper(), 1.7) for e in elements], dtype=np.float32)
    sphere = _fibonacci_sphere(n_points).astype(np.float32)
    pts = coords[:, None, :] + (radii[:, None, None] + probe) * sphere[None, :, :]  # (N,S,3)
    N = len(coords)
    tree = cKDTree(coords)
    pts_flat = pts.reshape(-1, 3)
    owner = np.repeat(np.arange(N), n_points)
    # a point is buried if any neighbor atom j (j != owner) has |p - c_j| < vdw_j + probe
    K = 16
    dist, idx = tree.query(pts_flat, k=K, workers=-1)
    idx = np.clip(idx, 0, N - 1)
    valid_nb = np.isfinite(dist)
    vdw_nb = radii[idx]                      # (P,K)
    owner_col = owner[:, None]
    buried = (dist < vdw_nb + probe) & (idx != owner_col) & valid_nb
    accessible = ~buried.any(1)
    frac = accessible.reshape(N, n_points).mean(1)
    # SASA convention: area of the probe-center surface (radius r_vdW + probe),
    # consistent with FreeSASA/Lee-Richards and the Tien et al. MAX_ASA reference values.
    sasa = 4 * np.pi * (radii + probe) ** 2 * frac
    return sasa, pts_flat, accessible

def run_dssp(pdb_path, dssp_exe="mkdssp"):
    """Parse classic DSSP output -> dict (chain, resnum) -> 3-state code."""
    out_path = pdb_path + ".dssp"
    subprocess.run([dssp_exe, "--output-format", "dssp", pdb_path, out_path],
                   capture_output=True, text=True, timeout=300)
    ss_map = {}
    started = False
    with open(out_path) as f:
        for line in f:
            if line.startswith("  #  RESIDUE"):
                started = True; continue
            if started and len(line) > 16:
                resnum = line[5:10].strip()
                if not resnum.isdigit():
                    continue
                ch = line[11]
                code = line[16]
                if code == " ":
                    code = "C"
                ss_map[(ch, int(resnum))] = {"H": "H", "B": "E", "E": "E"}.get(code, "C")
    os.unlink(out_path)
    return ss_map

def compute_static_features(cif_path, chains, dssp_exe="mkdssp"):
    """Return (list of residue dicts, per-residue CA coords). One dict per observed residue."""
    st = read_structure(cif_path)
    recs = chain_residues(st, chains)
    if len(recs) < 10:
        return [], None
    # flatten atoms of these chains in residue order
    coords, elements, atom_owner = [], [], []
    for ri, (chn, r) in enumerate(recs):
        for a in r:
            coords.append([a.pos.x, a.pos.y, a.pos.z])
            elements.append(a.element.name)
            atom_owner.append(ri)
    coords = np.asarray(coords, dtype=np.float32)
    atom_owner = np.asarray(atom_owner)
    sasa_atom, pts_flat, acc = shrake_rupley(coords, elements)
    res_sasa = np.bincount(atom_owner, weights=sasa_atom, minlength=len(recs))
    # depth: nearest accessible surface point per atom
    surf = pts_flat[acc]
    if len(surf) > 0:
        dmin, _ = cKDTree(surf).query(coords, k=1, workers=-1)
    else:
        dmin = np.zeros(len(coords))
    counts = np.bincount(atom_owner, minlength=len(recs))
    res_depth = np.bincount(atom_owner, weights=dmin, minlength=len(recs)) / np.maximum(counts, 1)
    # CA geometry
    ca = np.array([[r["CA"][0].pos.x, r["CA"][0].pos.y, r["CA"][0].pos.z] for _, r in recs])
    n = len(recs)
    tree = cKDTree(ca)
    pairs = tree.query_pairs(10.0, output_type="ndarray")
    wcn = np.zeros(n)
    contacts = np.zeros(n, dtype=int)
    if len(pairs):
        dd = np.linalg.norm(ca[pairs[:, 0]] - ca[pairs[:, 1]], axis=1)
        contrib = 1.0 / dd ** 2
        np.add.at(wcn, pairs[:, 0], contrib)
        np.add.at(wcn, pairs[:, 1], contrib)
        np.add.at(contacts, pairs[:, 0], 1)
        np.add.at(contacts, pairs[:, 1], 1)
    centroid = ca.mean(0)
    rad = np.linalg.norm(ca - centroid, axis=1)
    radial = rad / rad.max()
    L = n
    term = np.array([min(i, L - 1 - i) / max(L - 1, 1) for i in range(n)])
    rsa = np.array([res_sasa[i] / MAX_ASA.get(recs[i][1].name, 200) for i in range(n)])
    rsa = np.clip(rsa, 0, 1)
    # neighborhood means (10 A, excluding self)
    nb_rsa = np.zeros(n); nb_dep = np.zeros(n); nb_wcn = np.zeros(n)
    for i in range(n):
        nb = tree.query_ball_point(ca[i], 10.0)
        nb = [j for j in nb if j != i]
        if nb:
            nb_rsa[i] = rsa[nb].mean(); nb_dep[i] = res_depth[nb].mean(); nb_wcn[i] = wcn[nb].mean()
    # secondary structure via mkdssp (needs a PDB of the selected chains; build a
    # clean protein-only sub-structure with unique single-character chain names)
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    used = set()
    short_of = {}
    st_sub = gemmi.Structure()
    st_sub.spacegroup_hm = "P 1"
    sub_model = gemmi.Model("1")
    for chn in chains:
        nm = chn if len(chn) == 1 and chn not in used else None
        if nm is None:
            nm = next(a for a in alphabet if a not in used)
        used.add(nm)
        short_of[chn] = nm
        new_chain = gemmi.Chain(nm)
        for r in st[0][chn]:
            info = gemmi.find_tabulated_residue(r.name)
            if info is not None and info.is_amino_acid():
                new_chain.add_residue(r)
        sub_model.add_chain(new_chain)
    st_sub.add_model(sub_model)
    st_sub.setup_entities()
    pdb_path = None
    try:
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
            tmp.write(b"HEADER    PROTEIN\n" + st_sub.make_pdb_string().encode())
            pdb_path = tmp.name
        ss_map = run_dssp(pdb_path, dssp_exe)
        ss = [ss_map.get((short_of[chn], r.seqid.num), "C") for chn, r in recs]
    except Exception:
        ss = ["C"] * n
    finally:
        if pdb_path and os.path.exists(pdb_path):
            os.unlink(pdb_path)
    rows = []
    for i, (chn, r) in enumerate(recs):
        rows.append(dict(
            chain=chn, resnum=r.seqid.num, resname=r.name,
            aa=r.name if r.name in AA20 else "X",
            rsa=float(rsa[i]), depth=float(res_depth[i]), wcn=float(wcn[i]),
            ss=ss[i], radial=float(radial[i]), terminal=float(term[i]),
            nb_rsa=float(nb_rsa[i]), nb_depth=float(nb_dep[i]), nb_wcn=float(nb_wcn[i]),
            n_contacts=int(contacts[i]),
        ))
    return rows, ca
