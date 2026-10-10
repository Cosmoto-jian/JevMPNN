"""Frozen ligand classification filters (protocol section 6.1, inclusion_exclusion_rules.md).

First filter: AHoJ/BioLiP2 biological relevance (applied by AHoJ-DB itself).
Second filter: the fixed lists below + heavy-atom rule (>= 6 non-H atoms).
"""

# Single-atom metal / halide ions (E01)
SINGLE_ATOM_IONS = {
    "ZN","CA","MG","MN","FE","FE2","CU","CU1","HG","CD","CO","NI","NA","K","CL","BR","I","F",
    "LI","CS","BA","SR","AL","AU","AG","PB","SN","PT","PD","YB","TB","GD","SM","EU","DY","HO",
    "ER","TM","LU","RB","MO","W","V","CR","ZR","HF","TI","GA","IN","TL","BI","SB","BE","CE",
    "PR","ND","LA","SC","Y","NB","TA","RH","RU","OS","IR","CU3","EU3","CE3","ZR1","CUA","CU6",
    "CU8","CUA5","FEA","FEO","MER","MG8","MN3","MN5","MN6","NA6","NAO","OS4","PB2","PT4","SM4",
}

# Multi-atom inorganic / buffer / crystallization additives (E01)
ADDITIVES = {
    # inorganic anions/cations
    "SO4","PO4","NO3","CO3","OCN","NCO","CMO","AZI","IOD","BO4","WO4","SE4","SEO3","CRO","BCT",
    "O4B","SUCX","PER","OXY","O2","CO2","NO","H2S","NH4","N2H","H2O2","S2O3","SCN","CNO","CYN",
    # buffers
    "MES","ACT","ACE","CIT","TRS","EPE","PIP","PIPES","MOPS","CAPS","CHES","TES","BIS","ADA",
    "MAL","ACY","CAD","CAC","FMT","BU3","BTB","BES","TES","TMA","IMD","IMIDAZOLE","MPD","DMS",
    "DMF","DIO","EtOH","MEO","ACN","N8E","NH4AC",
    # cryoprotectants / polyols / PEG fragments
    "GOL","EDO","PEG","PGE","PG0","PG4","PG5","PG6","PG7","PG8","PE8","PEE","PEO","1PE","2PE",
    "P6G","P33","MRD","DTT","TCEP","BME","MET","ME0","ME1","2ET","ETX","PRO","IPA","BU1","BU2",
    "PAM","PRX","XPE","SPK","SMM","SME","SUC","TRE","LMT","L2B","OGA","OGP","NGP","C8E","C8N",
    "C7E","C6E","CHS","CHT","DDM","LDA","NMN" if False else "ZCM",
    # reducing agents / radicals
    "MSO","RDX","GLU" if False else "GLA",
}

# Nucleic acid residue ligands (DNA/RNA bases and common modified residues) (E01)
NUCLEIC = {
    "A","C","G","U","DA","DT","DG","DC","PSU","5MU","5MC","DHU","M2G","OMC","OMG","OMU","H2U",
    "C34","C31","U31","G7M","M1G","M1A","M7G","S4U","D5MU","BRU","5BU","5CU","5IU","BRG","2MG",
    "YYG","FHU","MOD","TH6","TPN","TP6","TCP","TDR","ADN","CYT","GUA","URA","THM","ARB","3TD",
}

# Amino acids and common modified amino acids (peptide/protein ligands) (E01)
AMINO_ACIDS = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO",
    "SER","THR","TRP","TYR","VAL","MSE","SEC","PYL","CSO","CME","CSS","OCS","HYP","MLY","MLZ",
    "M3L","ALY","KCX","LLP","SEP","TPO","PTR","THC","CAS","CSX","SNC","SMC","FME","ACE","NH2",
    "ABA","DAL","DLE","DNE","DVA","ORN","SAR","AIB","CGU","ASX","GLX","UNK","GLZ","TRQ","MEN",
}

def ligand_class(lig: str) -> str:
    """Return exclusion category or 'OK'."""
    L = lig.upper()
    if L in {"HOH", "DOD", "WAT"}:
        return "WATER"
    if L in SINGLE_ATOM_IONS:
        return "ION"
    if L in ADDITIVES:
        return "ADDITIVE"
    if L in NUCLEIC:
        return "NUCLEIC"
    if L in AMINO_ACIDS:
        return "PEPTIDE"
    return "OK"


def heavy_atoms_from_formula(formula: str) -> int:
    """Count non-hydrogen atoms from a CCD chemical formula like 'C21 H27 N7 O14 P2'."""
    import re
    total = 0
    for elem, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if not elem:
            continue
        if elem == "H":
            continue
        total += int(n) if n else 1
    return total
