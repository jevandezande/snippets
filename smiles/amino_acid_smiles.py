amino_acids = {
    "A": "C",
    "R": "CCCCNC(N)=N",
    "N": "CC(N)=O",
    "D": "CC(O)=O",
    "B": "CC(O)=O",
    "C": "CS",
    "E": "CCC(O)=O",
    "Q": "CCC(N)=O",
    "Z": "CCC(N)=O",
    "G": "[H]",
    "H": "CC1=CNC=N1",
    "I": "C(CC)([H])C",
    "L": "CC(C)C",
    "K": "CCCCN",
    "M": "CCSC",
    "F": "CC1=CC=CC=C1",
    "P": "C2CCCN2",
    "S": "CO",
    "T": "C(C)([H])O",
    "W": "CCC1=CNC2=C1C=CC=C2",
    "Y": "CC1=CC=C(O)C=C1",
    "V": "C(C)C",
}


full_smiles = {letter: f"OC(=O)C({smiles})N" for letter, smiles in amino_acids.items()}


def aa_to_smiles(aa):
    """
    >>> aa_to_smiles("A")
    'OC(=O)C(C)N'
    >>> aa_to_smiles("N")
    'OC(=O)C(CC(N)=O)N'
    """
    return full_smiles[aa]


def aas_to_smiles(aas):
    """
    >>> aas_to_smiles("CLB")
    'OC(=O)C(CS)NC(=O)C(CC(C)C)NC(=O)C(CC(O)=O)N'
    """
    return "O" + "".join(f"C(=O)C({amino_acids[aa]})N" for aa in aas)