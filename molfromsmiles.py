#!/usr/bin/env python3

import re
from decimal import Decimal
_masses_ = {
    "Ac": "227.0",
    "Al": "26.9815",
    "Am": "243.0",
    "Sb": "121.75",
    "Ar": "39.948",
    "As": "74.9216",
    "At": "210.0",
    "Ba": "137.0",
    "Bk": "247.0",
    "Be": "9.0122",
    "Bi": "208.98",
    "B": "10.81",
    "Br": "79.904",
    "Cd": "112.4",
    "Ca": "40.08",
    "Cf": "251.0",
    "C": "12.011",
    "Ce": "140.12",
    "Cs": "132.9054",
    "Cl": "35.453",
    "Cr": "51.996",
    "Co": "58.9332",
    "Cu": "63.546",
    "Cm": "247.0",
    "Dy": "162.5",
    "Es": "254.0",
    "Er": "167.26",
    "Eu": "151.96",
    "Fm": "257.0",
    "F": "18.9984",
    "Fr": "223.0",
    "Gd": "157.25",
    "Ga": "69.72",
    "Ge": "72.59",
    "Au": "196.966",
    "Hf": "178.49",
    "He": "4.0026",
    "Ho": "164.93",
    "H": "1.0079",
    "In": "114.82",
    "I": "126.904",
    "Ir": "192.22",
    "Fe": "55.847",
    "Kr": "83.8",
    "La": "138.905",
    "Lr": "256.0",
    "Pb": "207.2",
    "Li": "6.941",
    "Lu": "174.97",
    "Mg": "24.305",
    "Mn": "54.938",
    "Md": "258.0",
    "Hg": "200.59",
    "Mo": "95.94",
    "Nd": "144.24",
    "Ne": "20.179",
    "Np": "237.048",
    "Ni": "58.7",
    "Nb": "92.9064",
    "N": "14.0067",
    "No": "255.0",
    "Os": "190.2",
    "O": "15.9994",
    "Pd": "106.4",
    "P": "30.9738",
    "Pt": "195.09",
    "Pu": "244.0",
    "Po": "210.0",
    "K": "39.098",
    "Pr": "140.908",
    "Pm": "147.0",
    "Pa": "231.036",
    "Ra": "226.025",
    "Rn": "222.0",
    "Re": "186.207",
    "Rh": "102.906",
    "Rb": "85.4678",
    "Ru": "101.07",
    "Rf": "261.0",
    "Sm": "150.4",
    "Sc": "44.9559",
    "Se": "78.96",
    "Si": "28.086",
    "Ag": "107.868",
    "Na": "22.9898",
    "Sr": "87.62",
    "S": "32.06",
    "Ta": "180.948",
    "Tc": "98.9062",
    "Te": "127.6",
    "Tb": "158.925",
    "Tl": "204.37",
    "Th": "232.038",
    "Tm": "168.934",
    "Sn": "118.69",
    "Ti": "47.9",
    "W": "183.85",
    "U": "238.029",
    "V": "50.9414",
    "Xe": "131.3",
    "Yb": "173.04",
    "Y": "88.9059",
    "Zn": "65.38",
    "Zr": "91.22"
}
__pattern__ = "((?P<atom>({}))(?P<count>[1-9]*))".format("|".join(_masses_.keys()))


def smiles_to_mol(smiles_string):
    matches = re.finditer(__pattern__, smiles_string)
    atoms_and_counts = {}
    for match in matches:
        atom = match.group("atom")
        count = int(match.group("count")) if match.group("count") else 1
        if atom not in atoms_and_counts:
            atoms_and_counts[atom] = 0
        atoms_and_counts[atom] += count

    formula = ""
    total_mass = Decimal(0)
    for atom in sorted(atoms_and_counts.keys()):
        formula += atom
        count = atoms_and_counts[atom]
        if count > 1:
            formula += str(count)
        total_mass += count * Decimal(_masses_[atom])

    return {
        "formula": formula,
        "mass": total_mass
    }


if __name__ == "__main__":
    # calling this from the command line doesn't really makes sense,
    # since bash doesn't accept unescaped parentheses "(" ")" commonly
    # found in SMILES string
    from sys import argv
    mol = smiles_to_mol(argv[1])
    print("Formula: {} Atomic mass: {} g/mol".format(mol["formula"], mol["mass"]))
