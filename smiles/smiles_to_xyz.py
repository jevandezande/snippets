#!/usr/bin/env python

from argparse import ArgumentParser

from rdkit.Chem import AddHs, MolFromSmiles, MolToXYZFile
from rdkit.Chem.AllChem import EmbedMolecule, MMFFOptimizeMolecule


def smiles_to_xyz(smiles: str, output: str = "geom.xyz"):
    mol = MolFromSmiles(smiles)
    mol = AddHs(mol)

    EmbedMolecule(mol)
    MMFFOptimizeMolecule(mol)

    MolToXYZFile(mol, output)


def main():
    parser = ArgumentParser(description="")
    parser.add_argument("input", help="The smiles to be read", type=str)
    parser.add_argument("output", help="Where to output the xyz", type=str)

    args = parser.parse_args()

    smiles_to_xyz(args.input, args.output)


if __name__ == "__main__":
    main()
