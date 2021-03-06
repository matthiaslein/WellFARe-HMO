#!/usr/local/bin/python3

import argparse
import sys
import time
from importlib.util import find_spec

from src import wellfareSTO


def timestamp(s):
    print(s + time.strftime("%Y/%m/%d %X"))


# ASCII FONTS from: http://patorjk.com/software/taag/
# Font = "small"
def ProgramHeader():
    print("###########################################################")
    print("Wellington Fast Assessment of Reactions - Hückel MO Theory")
    print(" __      __   _ _ ___ _   ___     _  _ __  __  ___        ")
    print(" \ \    / /__| | | __/_\ | _ \___| || |  \/  |/ _ \       ")
    print("  \ \/\/ / -_) | | _/ _ \|   / -_) __ | |\/| | (_)|       ")
    print("   \_/\_/\___|_|_|_/_/ \_\_|_\___|_||_|_|  |_|\___/       ")
    print("                                              Version 0.01")
    print("    WellFAReHMO Copyright (C) 2015 Matthias Lein          ")
    print("   This program comes with ABSOLUTELY NO WARRANTY         ")
    print("    This is free software, and you are welcome to         ")
    print("      redistribute it under certain conditions.           ")
    timestamp('Program started at: ')
    print("###########################################################\n")


def ProgramFooter():
    print("\n###############################################################################")
    print("  ___                                ___         _    ")
    print(" | _ \_ _ ___  __ _ _ _ __ _ _ __   | __|_ _  __| |___")
    print(" |  _/ '_/ _ \/ _` | '_/ _` | '  \  | _|| ' \/ _` (_-<")
    print(" |_| |_| \___/\__, |_| \__,_|_|_|_| |___|_||_\__,_/__/")
    print("              |___/                                   ")
    timestamp('Program terminated at: ')
    print("###############################################################################")


def ProgramAbort():
    print("\n###############################################################################")
    print("  ___                  _             _          _ ")
    print(" | _ \_  _ _ _    __ _| |__  ___ _ _| |_ ___ __| |")
    print(" |   / || | ' \  / _` | '_ \/ _ \ '_|  _/ -_) _` |")
    print(" |_|_\\\_,_|_||_| \__,_|_.__/\___/_|  \__\___\__,_|")
    timestamp('Program aborted at: ')
    print("###############################################################################")
    sys.exit()
    return


def ProgramWarning(warntext=''):
    print("\n###############################################################################")
    print(" __      __             _           ")
    print(" \ \    / /_ _ _ _ _ _ (_)_ _  __ _ ")
    print("  \ \/\/ / _` | '_| ' \| | ' \/ _` |")
    print("   \_/\_/\__,_|_| |_||_|_|_||_\__, |")
    print("                              |___/ ")
    timestamp('Warning time/date: ')
    if warntext != '':
        print("###############################################################################")
        print("# ", warntext)
    print("###############################################################################\n")
    return


def ProgramError(errortext=''):
    print("\n###############################################################################")
    print("  ___                 ")
    print(" | __|_ _ _ _ ___ _ _ ")
    print(" | _|| '_| '_/ _ \ '_|")
    print(" |___|_| |_| \___/_|  ")
    timestamp('Error time/date: ')
    if errortext != '':
        print("###############################################################################")
        print("# ", errortext)
    print("###############################################################################")
    return


# Check for numpy, exit immediately if not available
module_loader = find_spec('numpy')
found = module_loader is not None
if not found:
    ProgramError("Module numpy is required")
    ProgramAbort()
import numpy as np

# Check for scipy, exit immediately if not available
module_loader = find_spec('scipy')
found = module_loader is not None
if not found:
    ProgramError("Module scipy is required")
    ProgramAbort()
import scipy.optimize


#############################################################################################################
# This section is for the definition of *all* constants and conversion factors
#############################################################################################################

# Conversion of mass in atomic mass units (AMU) to
# atomic units (electron masses)
def AMU2au(amu):
    return amu * 1822.88839


# Same in reverse
def au2AMU(au):
    return au / 1822.88839


# Conversion of length in Angstroms to  to
# atomic units (Bohrs)
def Ang2Bohr(ang):
    return ang * 1.889725989


# Same in reverse
def Bohr2Ang(bohr):
    return bohr / 1.889725989


SymbolToNumber = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9,
    "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17,
    "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24,
    "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31,
    "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38,
    "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45,
    "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51, "Te": 52,
    "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59,
    "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66,
    "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72, "Ta": 73,
    "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
    "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87,
    "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91, "U": 92, "Np": 93, "Pu": 94,
    "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100, "Md": 101,
    "No": 102, "Lr": 103, "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107,
    "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111, "Cn": 112, "Uut": 113,
    "Fl": 114, "Uup": 115, "Lv": 116, "Uus": 117, "Uuo": 118}

# Invert the above: atomic numbers to atomic symbols
NumberToSymbol = {v: k for k, v in SymbolToNumber.items()}

SymbolToMass = {
    "H": 1.00794, "He": 4.002602, "Li": 6.941, "Be": 9.012182, "B": 10.811,
    "C": 12.0107, "N": 14.0067, "O": 15.9994, "F": 18.9984032, "Ne": 20.1797,
    "Na": 22.98976928, "Mg": 24.3050, "Al": 26.9815386, "Si": 28.0855,
    "P": 30.973762, "S": 32.065, "Cl": 35.453, "Ar": 39.948, "K": 39.0983,
    "Ca": 40.078, "Sc": 44.955912, "Ti": 47.867, "V": 50.9415, "Cr": 51.9961,
    "Mn": 54.938045, "Fe": 55.845, "Co": 58.933195, "Ni": 58.6934, "Cu": 63.546,
    "Zn": 65.38, "Ga": 69.723, "Ge": 72.64, "As": 74.92160, "Se": 78.96,
    "Br": 79.904, "Kr": 83.798, "Rb": 85.4678, "Sr": 87.62, "Y": 88.90585,
    "Zr": 91.224, "Nb": 92.90638, "Mo": 95.96, "Tc": 98.0, "Ru": 101.07,
    "Rh": 102.90550, "Pd": 106.42, "Ag": 107.8682, "Cd": 112.411, "In": 114.818,
    "Sn": 118.710, "Sb": 121.760, "Te": 127.60, "I": 126.90447, "Xe": 131.293,
    "Cs": 132.9054519, "Ba": 137.327, "La": 138.90547, "Ce": 140.116,
    "Pr": 140.90765, "Nd": 144.242, "Pm": 145.0, "Sm": 150.36, "Eu": 151.964,
    "Gd": 157.25, "Tb": 158.92535, "Dy": 162.500, "Ho": 164.93032, "Er": 167.259,
    "Tm": 168.93421, "Yb": 173.054, "Lu": 174.9668, "Hf": 178.49, "Ta": 180.94788,
    "W": 183.84, "Re": 186.207, "Os": 190.23, "Ir": 192.217, "Pt": 195.084,
    "Au": 196.966569, "Hg": 200.59, "Tl": 204.3833, "Pb": 207.2, "Bi": 208.98040,
    "Po": 209.0, "At": 210.0, "Rn": 222.0, "Fr": 223.0, "Ra": 226.0, "Ac": 227.0,
    "Th": 232.03806, "Pa": 231.03588, "U": 238.02891, "Np": 237.0, "Pu": 244.0,
    "Am": 243.0, "Cm": 247.0, "Bk": 247.0, "Cf": 251.0, "Es": 252.0, "Fm": 257.0,
    "Md": 258.0, "No": 259.0, "Lr": 262.0, "Rf": 267.0, "Db": 268.0, "Sg": 271.0,
    "Bh": 272.0, "Hs": 270.0, "Mt": 276.0, "Ds": 281.0, "Rg": 280.0, "Cn": 285.0,
    "Uut": 284.0, "Uuq": 289.0, "Uup": 288.0, "Uuh": 293.0, "Uuo": 294.0}

# Define dictionary to convert atomic symbols to covalent radii (in Angstrom)
SymbolToRadius = {
    "H": 0.37, "He": 0.32, "Li": 1.34, "Be": 0.90, "B": 0.82, "C": 0.77,
    "N": 0.75, "O": 0.73, "F": 0.71, "Ne": 0.69, "Na": 1.54, "Mg": 1.30,
    "Al": 1.18, "Si": 1.11, "P": 1.06, "S": 1.02, "Cl": 0.99, "Ar": 0.97,
    "K": 1.96, "Ca": 1.74, "Sc": 1.44, "Ti": 1.36, "V": 1.25, "Cr": 1.27,
    "Mn": 1.39, "Fe": 1.25, "Co": 1.26, "Ni": 1.21, "Cu": 1.38, "Zn": 1.31,
    "Ga": 1.26, "Ge": 1.22, "As": 1.19, "Se": 1.16, "Br": 1.14, "Kr": 1.10,
    "Rb": 2.11, "Sr": 1.92, "Y": 1.62, "Zr": 1.48, "Nb": 1.37, "Mo": 1.45,
    "Tc": 1.56, "Ru": 1.26, "Rh": 1.35, "Pd": 1.31, "Ag": 1.53, "Cd": 1.48,
    "In": 1.44, "Sn": 1.41, "Sb": 1.38, "Te": 1.35, "I": 1.33, "Xe": 1.30,
    "Cs": 2.25, "Ba": 1.98, "La": 1.69, "Ce": 1.70, "Pr": 1.70, "Nd": 1.70,
    "Pm": 1.70, "Sm": 1.70, "Eu": 1.70, "Gd": 1.70, "Tb": 1.70, "Dy": 1.70,
    "Ho": 1.70, "Er": 1.70, "Tm": 1.70, "Yb": 1.70, "Lu": 1.60, "Hf": 1.50,
    "Ta": 1.38, "W": 1.46, "Re": 1.59, "Os": 1.28, "Ir": 1.37, "Pt": 1.28,
    "Au": 1.44, "Hg": 1.49, "Tl": 1.48, "Pb": 1.47, "Bi": 1.46, "Po": 1.50,
    "At": 1.50, "Rn": 1.45, "Fr": 1.50, "Ra": 1.50, "Ac": 1.50, "Th": 1.50,
    "Pa": 1.50, "U": 1.50, "Np": 1.50, "Pu": 1.50, "Am": 1.50, "Cm": 1.50,
    "Bk": 1.50, "Cf": 1.50, "Es": 1.50, "Fm": 1.50, "Md": 1.50, "No": 1.50,
    "Lr": 1.50, "Rf": 1.50, "Db": 1.50, "Sg": 1.50, "Bh": 1.50, "Hs": 1.50,
    "Mt": 1.50, "Ds": 1.50, "Rg": 1.50, "Cn": 1.50, "Uut": 1.50, "Uuq": 1.50,
    "Uup": 1.50, "Uuh": 1.50, "Uus": 1.50, "Uuo": 1.50}

# Define dictionary to convert atomic symbols to (Pauling) electronegativity
SymbolToEN = {
    "H": 2.20, "He": 0.00, "Li": 0.98, "Be": 1.57, "B": 2.04, "C": 2.55,
    "N": 3.04, "O": 3.44, "F": 3.98, "Ne": 0.00, "Na": 0.93, "Mg": 1.31,
    "Al": 1.61, "Si": 1.90, "P": 2.19, "S": 2.58, "Cl": 3.16, "Ar": 0.00,
    "K": 0.82, "Ca": 1.00, "Sc": 1.36, "Ti": 1.54, "V": 1.63, "Cr": 1.66,
    "Mn": 1.55, "Fe": 1.83, "Co": 1.88, "Ni": 1.91, "Cu": 1.90, "Zn": 1.65,
    "Ga": 1.81, "Ge": 2.01, "As": 2.18, "Se": 2.55, "Br": 2.96, "Kr": 3.00,
    "Rb": 0.82, "Sr": 0.95, "Y": 1.22, "Zr": 1.33, "Nb": 1.60, "Mo": 2.16,
    "Tc": 1.90, "Ru": 2.00, "Rh": 2.28, "Pd": 2.20, "Ag": 1.93, "Cd": 1.69,
    "In": 1.78, "Sn": 1.96, "Sb": 2.05, "Te": 2.10, "I": 2.66, "Xe": 2.60,
    "Cs": 0.79, "Ba": 0.89, "La": 1.10, "Ce": 1.12, "Pr": 1.13, "Nd": 1.14,
    "Pm": 1.13, "Sm": 1.17, "Eu": 1.20, "Gd": 1.20, "Tb": 1.10, "Dy": 1.22,
    "Ho": 1.23, "Er": 1.24, "Tm": 1.25, "Yb": 1.10, "Lu": 1.27, "Hf": 1.30,
    "Ta": 1.50, "W": 2.36, "Re": 1.90, "Os": 2.20, "Ir": 2.20, "Pt": 2.28,
    "Au": 2.54, "Hg": 2.00, "Tl": 1.62, "Pb": 1.87, "Bi": 2.02, "Po": 2.00,
    "At": 2.20, "Rn": 2.20, "Fr": 0.70, "Ra": 0.90, "Ac": 1.10, "Th": 1.30,
    "Pa": 1.50, "U": 1.38, "Np": 1.36, "Pu": 1.28, "Am": 1.13, "Cm": 1.28,
    "Bk": 1.30, "Cf": 1.30, "Es": 1.30, "Fm": 1.30, "Md": 1.30, "No": 1.30,
    "Lr": 1.30, "Rf": 1.30, "Db": 1.30, "Sg": 1.30, "Bh": 1.30, "Hs": 1.30,
    "Mt": 1.30, "Ds": 1.30, "Rg": 1.30, "Cn": 1.30, "Uut": 1.30, "Uuq": 1.30,
    "Uup": 1.30, "Uuh": 1.30, "Uus": 1.30, "Uuo": 1.30}

# Define dictionary to convert atomic symbols to number of valence electrons
SymbolToValE = {
    "H": 1, "He": 2, "Li": 1, "Be": 2, "B": 3, "C": 4, "N": 5, "O": 6, "F": 7, "Ne": 8,
    "Na": 1, "Mg": 2, "Al": 3, "Si": 4, "P": 5, "S": 6, "Cl": 7, "Ar": 8
}

# Define dictionary to convert atomic symbols to STO exponents for s functions
SymbolToSTOexpS = {
    "H": 1.200, "He": 1.688, "Li": 0.650, "Be": 0.975, "B": 1.300, "C": 1.625,
    "N": 1.950, "O": 2.275, "F": 2.425, "Ne": 2.879, "Na": 0.733, "Mg": 1.100,
    "Al": 1.167, "Si": 1.383, "P": 1.750, "S": 2.122, "Cl": 2.183, "Ar": 2.461,

}

# Define dictionary to convert atomic symbols to STO exponents for p functions
SymbolToSTOexpP = {
    "Li": 0.650, "Be": 0.975, "B": 1.300, "C": 1.625,
    "N": 1.950, "O": 2.275, "F": 2.425, "Ne": 2.879, "Na": 0.733, "Mg": 1.100,
    "Al": 1.167, "Si": 1.383, "P": 1.300, "S": 1.827, "Cl": 1.733, "Ar": 2.105,
}

# Define dictionary to convert atomic symbols to ionisation energies in Extended Hückel Hamiltonians
# Here for s electrons in hartrees
SymbolToEHTieS = {
    "H": -0.5000, "He": -0.8599, "Li": -0.1984, "Be": -0.3675, "B": -0.5586, "C": -0.7144,
    "N": -0.9555, "O": -1.1870, "F": -1.4700, "Ne": -1.5876, "Na": -0.1874, "Mg": -0.3307,
    "Al": -0.4520, "Si": -0.6358, "P": -0.6835, "S": -0.7350, "Cl": -0.9665, "Ar": -1.2774,
}

# Define dictionary to convert atomic symbols to ionisation energies in Extended Hückel Hamiltonians
# Here for p electrons in hartrees
SymbolToEHTieP = {
    "Li": -0.1286, "Be": -0.2205, "B": -0.3124, "C": -0.3921,
    "N": -0.4924, "O": -0.5439, "F": -0.6652, "Ne": -0.7350, "Na": -0.1102, "Mg": -0.1654,
    "Al": -0.2389, "Si": -0.3381, "P": -0.5145, "S": -0.4042, "Cl": -0.5218, "Ar": -0.5910,
}

# Define dictionary to convert angular momentum quantum numbers to symbols
L2Symb = {
    0: "s", 1: "p", 2: "d", 3: "f", 4: "g", 5: "h", 6: "i", 7: "j", 8: "k"
}

# Define dictionary to convert quantum numbers to symbols for p orbitals
Porb2Symb = {
    -1: "x", 0: "z", 1: "y"
}

# Define dictionary to convert quantum numbers to symbols for d orbitals
Dorb2Symb = {
    -2: "xx", -1: "xx", 0: "z2", 1: "yy", 2: "yy"
}


# Function to translate quantum numbers to symbols.
def qn2symb(l, m=None):
    if m == None:
        return L2Symb[l]
    elif l == 0:
        return ""
    elif l == 1:
        return Porb2Symb[m]
    elif l == 2:
        return Dorb2Symb[m]


#############################################################################################################
# Do *not* define constants or conversion factors below here
#############################################################################################################

# Test if the argument is (can be converted to)
# an integer number
def isInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


#############################################################################################################
# STO (Slater Type Orbital class and class methods to be defined below
#############################################################################################################

class STO:
    """ A Slater Type Orbital with an atomic symbol quantum numbers n and l, and an exponent"""

    def __init__(self, sym, n, l, exp=None, ie=None):
        """ (STO, str, number, number, number, number) -> NoneType

        Create an STO with (int) quantum numbers n and l.
        Exponent exp and ionisation energy are set automatically according to
        symbol if not explicitly specified.
        """

        self.n = n
        self.l = l
        if l == 0:
            self.exp = SymbolToSTOexpS[sym]
            self.ie = SymbolToEHTieS[sym]
        elif l == 1:
            self.exp = SymbolToSTOexpP[sym]
            self.ie = SymbolToEHTieP[sym]
        else:
            ProgramError("This angular momentum is not (yet) implemented")
            ProgramAbort()

    def __str__(self):
        """ (STO) -> str

        Return a string representation of this STO in this format:

          (n, l, exp, ie)
        """

        return '({0}, {1}, {2}, {3})'.format(self.n, self.l, self.exp, self.ie)

    def __repr__(self):
        """ (STO) -> str

        Return a string representation of this STO in this format:"

          STO(n, l, exp, ie)
        """

        return '({0}, {1}, {2}, {3})'.format(self.n, self.l, self.exp, self.ie)


#############################################################################################################
# Atom class and class methods to be defined below
#############################################################################################################

class Atom:
    """ An atom with an atomic symbol and cartesian coordinates"""

    def __init__(self, sym, x, y, z, basis="eht"):
        """ (Atom, int, str, number, number, number) -> NoneType

        Create an Atom with (string) symbol sym,
        and (float) cartesian coordinates (x, y, z).
        (int) charge charge and (float) mass mass are set
        automatically according to symbol. STO basis functions
        are added automatically according to choice
        """

        self.symbol = sym
        self.charge = SymbolToNumber[sym]
        self.mass = SymbolToMass[sym]
        self.coord = [x, y, z]
        self.basis = []
        self.valele = SymbolToValE[sym]

        if sym == "H":
            self.basis.append(STO(sym, 1, 0))
        elif sym == "He":
            self.basis.append(STO(sym, 1, 0))
        elif sym == "Li":
            self.basis.append(STO(sym, 2, 0))
            self.basis.append(STO(sym, 2, 1))
        elif sym == "Be":
            self.basis.append(STO(sym, 2, 0))
            self.basis.append(STO(sym, 2, 1))
        elif sym == "B":
            self.basis.append(STO(sym, 2, 0))
            self.basis.append(STO(sym, 2, 1))
        elif sym == "C":
            self.basis.append(STO(sym, 2, 0))
            self.basis.append(STO(sym, 2, 1))
        elif sym == "N":
            self.basis.append(STO(sym, 2, 0))
            self.basis.append(STO(sym, 2, 1))
        elif sym == "O":
            self.basis.append(STO(sym, 2, 0))
            self.basis.append(STO(sym, 2, 1))
        elif sym == "F":
            self.basis.append(STO(sym, 2, 0))
            self.basis.append(STO(sym, 2, 1))
        elif sym == "Ne":
            self.basis.append(STO(sym, 2, 0))
            self.basis.append(STO(sym, 2, 1))
        elif sym == "Na":
            self.basis.append(STO(sym, 3, 0))
            self.basis.append(STO(sym, 3, 1))
        elif sym == "Mg":
            self.basis.append(STO(sym, 3, 0))
            self.basis.append(STO(sym, 3, 1))
        elif sym == "Al":
            self.basis.append(STO(sym, 3, 0))
            self.basis.append(STO(sym, 3, 1))
        elif sym == "Si":
            self.basis.append(STO(sym, 3, 0))
            self.basis.append(STO(sym, 3, 1))
        elif sym == "P":
            self.basis.append(STO(sym, 3, 0))
            self.basis.append(STO(sym, 3, 1))
        elif sym == "S":
            self.basis.append(STO(sym, 3, 0))
            self.basis.append(STO(sym, 3, 1))
        elif sym == "Cl":
            self.basis.append(STO(sym, 3, 0))
            self.basis.append(STO(sym, 3, 1))
        elif sym == "Ar":
            self.basis.append(STO(sym, 3, 0))
            self.basis.append(STO(sym, 3, 1))

    def __str__(self):
        """ (Atom) -> str

        Return a string representation of this Atom in this format:

          (SYMBOL, X, Y, Z, [basis])
        """
        s = ""
        for i in self.basis:
            s = s + i.__str__()
            if i != self.basis[-1]:
                s = s + " "

        return '({0}, ({1}, {2}, {3}), ({4})'.format(self.symbol, self.coord[0], self.coord[1], self.coord[2], s)

    def __repr__(self):
        """ (Atom) -> str

        Return a string representation of this Atom in this format:"

          Atom("SYMBOL", charge, mass, X, Y, Z)
        """

        return '("{0}", {1}, {2}, {3}, {4}, {5})'.format(self.symbol, self.charge, self.mass, self.coord[0],
                                                         self.coord[1], self.coord[2])

    def x(self, x):
        """ (Atom) -> NoneType

        Set x coordinate to x
        """
        self.coord[0] = x

    def y(self, y):
        """ (Atom) -> NoneType

        Set y coordinate to y
        """
        self.coord[1] = y

    def z(self, z):
        """ (Atom) -> NoneType

        Set z coordinate to z
        """
        self.coord[2] = z


#############################################################################################################
# Molecule class and class methods to be defined below
#############################################################################################################

class Molecule:
    """A molecule with a name, charge and a list of atoms"""

    def __init__(self, name, charge=0):
        """ (Molecule, str, int) -> NoneType

        Create a Molecule named name with charge charge and no atoms
        (int) Multiplicity mult is automatically set to the lowest
        possible value (1 or 2).
        """

        self.name = name
        self.charge = charge
        self.mult = 1
        self.atoms = []

    def addAtom(self, a, verbosity=0):
        """ (Molecule, Atom) -> NoneType

        Add a to my list of Atoms.
        """

        self.atoms.append(a)
        if verbosity >= 2:
            print(" {:<3} {: .8f} {: .8f} {: .8f}".format(a.symbol, a.coord[0], a.coord[1], a.coord[2]))
        nucchg = 0
        for i in self.atoms:
            nucchg = nucchg + i.charge
        if (nucchg - self.charge) % 2 != 0:
            self.mult = 2
        else:
            self.mult = 1

    def __str__(self):
        """ (Molecule) -> str

        Return a string representation of this Molecule in this format:
        (NAME, CHARGE, MULT, (ATOM1, ATOM2, ...))
        """

        res = ''
        for atom in self.atoms:
            res = res + str(atom) + ', '
        res = res[:-2]
        return '({0}, {1}, {2}, ({3}))'.format(self.name, self.charge, self.mult, res)

    def __repr__(self):
        """ (Molecule) -> str

        Return a string representation of this Molecule in this format:
        (NAME, CHARGE, MULT, (ATOM1, ATOM2, ...))
        """

        res = ''
        for atom in self.atoms:
            res = res + str(atom) + ', '
        res = res[:-2]
        return '({0}, {1}, {2}, ({3}))'.format(self.name, self.charge, self.mult, res)

    def mass(self):
        """ (Molecule) -> number

        Return the molar mass as sum of atomic masses
        """

        mass = 0.0
        for atom in self.atoms:
            mass = mass + atom.mass

        return mass

    def numatoms(self):
        """ (Molecule) -> int

        Return the number of atoms in the molecule
        """

        return int(len(self.atoms))

    def chatom(self, n, at):
        """ (Molecule) -> NoneType

        Change the nth atom of the Molecule
        """

        self.atoms[n] = at

    def movatom(self, n, x, y, z):
        """ (Molecule) -> NoneType

        Move the nth atom to position x, y, z
        """

        self.atoms[n].x(x)
        self.atoms[n].y(y)
        self.atoms[n].z(z)

    def atmmass(self, n, m):
        """ (Molecule) -> NoneType

        Change the mass of the nth atom to m
        """

        self.atoms[n].mass = m

    def orient(self, verbosity=0):
        """ (Molecule) -> NoneType

        Translate centre of mass to coordinate origin and
        (re-)Orient the molecule along the principal axes of inertia.
        """

        # The molecular center of mass
        xValue = 0.0
        yValue = 0.0
        zValue = 0.0
        for i in self.atoms:
            xValue = xValue + (i.mass * i.coord[0])
            yValue = yValue + (i.mass * i.coord[1])
            zValue = zValue + (i.mass * i.coord[2])
        xValue = xValue / (self.mass())
        yValue = yValue / (self.mass())
        zValue = zValue / (self.mass())

        # Translate whole molecule into the center of mass reference frame
        for i in self.atoms:
            i.coord[0] = i.coord[0] - xValue
            i.coord[1] = i.coord[1] - yValue
            i.coord[2] = i.coord[2] - zValue

        # Build inertia tensor
        inertiaTensor = []
        Ixx = 0.0
        Ixy = 0.0
        Ixz = 0.0
        Iyx = 0.0
        Iyy = 0.0
        Iyz = 0.0
        Izx = 0.0
        Izy = 0.0
        Izz = 0.0
        for i in self.atoms:
            Ixx = Ixx + (
                i.mass * (
                    (Ang2Bohr(i.coord[1]) * Ang2Bohr(i.coord[1])) + (Ang2Bohr(i.coord[2]) * Ang2Bohr(i.coord[2]))))
            Ixy = Ixy - i.mass * Ang2Bohr(i.coord[0]) * Ang2Bohr(i.coord[1])
            Ixz = Ixz - i.mass * Ang2Bohr(i.coord[0]) * Ang2Bohr(i.coord[2])
            Iyx = Iyx - i.mass * Ang2Bohr(i.coord[1]) * Ang2Bohr(i.coord[0])
            Iyy = Iyy + (
                i.mass * (
                    (Ang2Bohr(i.coord[0]) * Ang2Bohr(i.coord[0])) + (Ang2Bohr(i.coord[2]) * Ang2Bohr(i.coord[2]))))
            Iyz = Iyz - i.mass * Ang2Bohr(i.coord[1]) * Ang2Bohr(i.coord[2])
            Izx = Izx - i.mass * Ang2Bohr(i.coord[2]) * Ang2Bohr(i.coord[0])
            Izy = Izy - i.mass * Ang2Bohr(i.coord[2]) * Ang2Bohr(i.coord[1])
            Izz = Izz + (
                i.mass * (
                    (Ang2Bohr(i.coord[0]) * Ang2Bohr(i.coord[0])) + (Ang2Bohr(i.coord[1]) * Ang2Bohr(i.coord[1]))))
        inertiaTensor.append([Ixx, Ixy, Ixz])
        inertiaTensor.append([Iyx, Iyy, Iyz])
        inertiaTensor.append([Izx, Izy, Izz])
        inertiaTensor = np.matrix(inertiaTensor)

        # Diagonalise inertia tensor
        inertiaMoments, inertialAxes = np.linalg.eig(inertiaTensor)

        # Orthogonalise eigenvectors (only sometimes necessary)...
        inertialAxes, r = np.linalg.qr(inertialAxes)

        # Sort moments from highest to lowest
        idx = inertiaMoments.argsort()[::-1]
        inertiaMoments = inertiaMoments[idx]
        inertialAxes = inertialAxes[:, idx]

        # Transform molecular coordinates into new frame of principal axes of inertia
        for i in self.atoms:
            vector = [i.coord[0], i.coord[1], i.coord[2]]
            vector = np.matrix(vector)
            vector = np.matrix.transpose(inertialAxes).dot(np.matrix.transpose(vector))
            vector = np.array(vector).flatten().tolist()
            i.coord[0] = vector[0]
            i.coord[1] = vector[1]
            i.coord[2] = vector[2]

        if verbosity == 1:
            print("\nRe-orienting molecule: {}".format(self.name))
        elif verbosity == 2:
            print("\nRe-orienting molecule: {}:".format(self.name))
            print("Moving center of mass to coordinate origin and aligning principal axes with coordinate axes.")
        elif verbosity == 3:
            print("\nRe-orienting molecule: {}:".format(self.name))
            print("Moving center of mass to coordinate origin and aligning principal axes with coordinate axes.")
            print("New coordinates:")
            print(self.gaussString())

    def cartesianCoordinates(self):
        """ (Molecule) ->

        Returns a list containing the cartesian coordinates.
        """

        coord = []
        for i in self.atoms:
            coord.append(i.coord[0])
            coord.append(i.coord[1])
            coord.append(i.coord[2])

        return coord

    def xyzString(self):
        """ (Molecule) -> str

        Returns a string containing an xyz file
        """

        s = str(self.numatoms()) + "\n" + self.name + "\n"
        for i in self.atoms:
            t = "{:<3} {: .8f} {: .8f} {: .8f}\n".format(i.symbol, i.coord[0], i.coord[1], i.coord[2])
            s = s + t

        return s

    def gamessString(self):
        """ (Molecule) -> str

          Returns a string containing cartesian coordinates in Gamess format
        """

        s = " $DATA\n" + self.name + "\nC1\n"
        for i in self.atoms:
            t = "{:<3} {:<3d} {: .8f} {: .8f} {: .8f}\n".format(i.symbol, i.charge, i.coord[0], i.coord[1], i.coord[2])
            s = s + t
        s = s + " $END\n"
        return s

    def gaussString(self):
        """ (Molecule) -> str

          Returns a string containing cartesian coordinates in Gaussian format
        """

        s = "\n" + self.name + "\n\n" + str(self.charge) + " " + str(self.mult) + "\n"
        for i in self.atoms:
            t = "{:<3} {: .8f} {: .8f} {: .8f}\n".format(i.symbol, i.coord[0], i.coord[1], i.coord[2])
            s = s + t
        s = s + "\n"
        return s

    def HMOEnergy(self, cartCoordinates, K=1.75, charge=0, verbosity=0):
        """ (Molecule) -> number (extended Hueckel aka Tight Binding energy)

          Returns a number containing the molecular energy according to the current extended Hueckel aka Tight Binding
          definition at structure specified by the provided cartesian coordinates.
        """

        # Assemble an array that holds information about the basis set.
        molbasis = []
        valence_electrons = 0
        for atomnum, i in enumerate(self.atoms):
            valence_electrons += i.valele
            for j in i.basis:
                for k in range(-1 * j.l, j.l + 1):
                    molbasis.append([atomnum, j.n, j.l, k, j.exp, j.ie])

        # Print the atomic basis of the calculation in a pretty way (with symbols instead of pure quantum numbers
        if verbosity >= 2 and verbosity < 3:
            print("\nBasis Functions")
            print("  Atom            exp     VSIE")
            for i in molbasis:
                print(
                    " {: >3}({: >3}) {:>2}{}{:<2}  {:.4f}  {: .5f}".format(self.atoms[i[0]].symbol, i[0], i[1],
                                                                           qn2symb(i[2]), qn2symb(i[2], i[3]), i[4],
                                                                           i[5]))
        # Print the atomic basis of the calculation with the quantum numbers themselves shown
        elif verbosity >= 3:
            print("\nBasis Functions")
            print("  Atom      n   l   m   exp     VSIE")
            for i in molbasis:
                print(
                    " {: >3}({: >3}) {: >3} {: >3} {: >3}  {:.4f}  {: .5f}".format(self.atoms[i[0]].symbol, i[0], i[1],
                                                                                   i[2], i[3], i[4], i[5]))

        # Create overlap matrix
        overlap = np.zeros((len(molbasis), len(molbasis)))
        # Calculate overlap matrix elements
        for i in range(0, len(molbasis)):
            # Exploit matrix symmetry by only calculating diagonal and upper triangle, then copying elements
            # to fill the rest
            for j in range(i, len(molbasis)):
                overlap[i][j] = wellfareSTO.SlaterOverlapCartesian(molbasis[i][1], molbasis[i][2], molbasis[i][3],
                                                                   molbasis[i][4],
                                                                   self.atoms[molbasis[i][0]].coord[0],
                                                                   self.atoms[molbasis[i][0]].coord[1],
                                                                   self.atoms[molbasis[i][0]].coord[2],
                                                                   molbasis[j][1],
                                                                   molbasis[j][2], molbasis[j][3], molbasis[j][4],
                                                                   self.atoms[molbasis[j][0]].coord[0],
                                                                   self.atoms[molbasis[j][0]].coord[1],
                                                                   self.atoms[molbasis[j][0]].coord[2])
                overlap[j][i] = overlap[i][j]
        if verbosity >= 3:
            # Print routine for the overlap matrix
            print("\nOverlap Matrix")
            s = [""] * (len(overlap) + 1)
            for i in range(0, len(overlap)):
                s[0] += "{: >3}({: >3}){:>2}{}{:<2}".format(self.atoms[molbasis[i][0]].symbol,
                                                            molbasis[i][0], molbasis[i][1],
                                                            qn2symb(molbasis[i][2]),
                                                            qn2symb(molbasis[i][2], molbasis[i][3]))
                for j in range(0, len(overlap)):
                    s[i + 1] += "   {: .6f} ".format(overlap[i][j])
            for i in range(0, len(s[0]), 65):
                for j in range(0, (len(overlap) + 1)):
                    if len(s[j]) < (i + 65):
                        if j == 0:
                            print("              " + s[j][i:len(s[0])])

                        else:
                            print("{: >3}({: >3}){:>2}{}{:<2}".format(self.atoms[molbasis[j - 1][0]].symbol,
                                                                      molbasis[j - 1][0], molbasis[j - 1][1],
                                                                      qn2symb(molbasis[j - 1][2]),
                                                                      qn2symb(molbasis[j - 1][2], molbasis[j - 1][3])) +
                                  s[j][i:len(s[0])])

                    else:
                        if j == 0:
                            print("              " + s[j][i:i + 65])
                        else:
                            print("{: >3}({: >3}){:>2}{}{:<2}".format(self.atoms[molbasis[j - 1][0]].symbol,
                                                                      molbasis[j - 1][0], molbasis[j - 1][1],
                                                                      qn2symb(molbasis[j - 1][2]),
                                                                      qn2symb(molbasis[j - 1][2], molbasis[j - 1][3])) +
                                  s[j][i:i + 65])
                print("")


        # Create Hamiltonian matrix
        hamiltonian = np.zeros((len(molbasis), len(molbasis)))
        # Calculate Hamiltonian matrix elements
        for i in range(0, len(molbasis)):
            # Exploit matrix symmetry by only calculating diagonal and upper triangle, then copying elements
            # to fill the rest
            for j in range(i, len(molbasis)):
                if i == j:
                    # Use Valence State Ionisation Energies for diagonal elements
                    hamiltonian[i][j] = molbasis[i][5]
                else:
                    # Use Wolfsberg-Helmholtz for off-diagonal elements
                    hamiltonian[i][j] = K * overlap[i][j] * ((molbasis[i][5] + molbasis[j][5]) / 2)
                    hamiltonian[j][i] = hamiltonian[i][j]
        if verbosity >= 3:
            # Print routine for the Hamiltonian matrix
            print("\nHamiltonian Matrix")
            s = [""] * (len(hamiltonian) + 1)
            for i in range(0, len(hamiltonian)):
                s[0] += "{: >3}({: >3}){:>2}{}{:<2}".format(self.atoms[molbasis[i][0]].symbol,
                                                            molbasis[i][0], molbasis[i][1],
                                                            qn2symb(molbasis[i][2]),
                                                            qn2symb(molbasis[i][2], molbasis[i][3]))
                for j in range(0, len(hamiltonian)):
                    s[i + 1] += "   {: .6f} ".format(hamiltonian[i][j])
            for i in range(0, len(s[0]), 65):
                for j in range(0, (len(hamiltonian) + 1)):
                    if len(s[j]) < (i + 65):
                        if j == 0:
                            print("              " + s[j][i:len(s[0])])

                        else:
                            print("{: >3}({: >3}){:>2}{}{:<2}".format(self.atoms[molbasis[j - 1][0]].symbol,
                                                                      molbasis[j - 1][0], molbasis[j - 1][1],
                                                                      qn2symb(molbasis[j - 1][2]),
                                                                      qn2symb(molbasis[j - 1][2], molbasis[j - 1][3])) +
                                  s[j][i:len(s[0])])

                    else:
                        if j == 0:
                            print("              " + s[j][i:i + 65])
                        else:
                            print("{: >3}({: >3}){:>2}{}{:<2}".format(self.atoms[molbasis[j - 1][0]].symbol,
                                                                      molbasis[j - 1][0], molbasis[j - 1][1],
                                                                      qn2symb(molbasis[j - 1][2]),
                                                                      qn2symb(molbasis[j - 1][2], molbasis[j - 1][3])) +
                                  s[j][i:i + 65])
                print("")

        # Use SciPy algorithm for generalised eigenvalue problem for symmetric matrices to solve
        # HC = SCE, H and S are our input matrices, E holds the energies and C are the coefficients.
        MOEnergies, MOVectors = scipy.linalg.eigh(hamiltonian, b=overlap)

        # Calculate total energy as sum over energies of occupied MOs
        energy = 0.0
        for i in range(0, valence_electrons):
            energy += MOEnergies[i // 2]

        # Print MO energies
        if verbosity >= 3:
            print("\nMO Energies ({} electrons, total energy {: .5f} hartree)".format(valence_electrons, energy))
            s = ""
            for i in range(0, len(MOEnergies)):
                s += " {: .5f}".format(MOEnergies[i])
            for i in range(0, len(s), 72):
                if len(s) < (i + 72):
                    print(s[i:len(s)])
                else:
                    print(s[i:i + 72])

        # Print MO vectors in a pretty way with MO numbers, energies and occupations. Each row is prefixed by the
        # atomic orbital that controbutes to the MO.
        if verbosity >= 2:
            print("\nMO Vectors")
            s = [""] * (len(MOVectors) + 3)
            for i in range(0, len(MOVectors)):
                s[0] += " {:>8} ".format(i + 1)
                s[1] += " {: .6f}".format(MOEnergies[i])
                if (i * 2) + 1 < valence_electrons:
                    s[2] += " {:>8} ".format(2)
                elif (i * 2) + 1 == valence_electrons:
                    s[2] += " {:>8} ".format(1)
                else:
                    s[2] += " {:>8} ".format(0)
                for j in range(0, len(MOVectors)):
                    s[i + 3] += " {: .6f}".format(MOVectors[i][j])
            for i in range(0, len(s[0]), 60):
                for j in range(0, (len(MOVectors) + 3)):
                    if len(s[j]) < (i + 60):
                        if j == 0:
                            print(" MO number:  " + s[j][i:len(s[0])])
                        elif j == 1:
                            print(" MO energy:  " + s[j][i:len(s[0])])
                        elif j == 2:
                            print(" MO occ   :  " + s[j][i:len(s[0])])
                        else:
                            print("{: >3}({: >3}){:>2}{}{:<2}".format(self.atoms[molbasis[j - 3][0]].symbol,
                                                                      molbasis[j - 3][0], molbasis[j - 3][1],
                                                                      qn2symb(molbasis[j - 3][2]),
                                                                      qn2symb(molbasis[j - 3][2], molbasis[j - 3][3])) +
                                  s[j][i:len(s[0])])
                    else:
                        if j == 0:
                            print(" MO number:  " + s[j][i:i + 60])
                        elif j == 1:
                            print(" MO energy:  " + s[j][i:i + 60])
                        elif j == 2:
                            print(" MO occ   :  " + s[j][i:i + 60])
                        else:
                            print("{: >3}({: >3}){:>2}{}{:<2}".format(self.atoms[molbasis[j - 3][0]].symbol,
                                                                      molbasis[j - 3][0], molbasis[j - 3][1],
                                                                      qn2symb(molbasis[j - 3][2]),
                                                                      qn2symb(molbasis[j - 3][2], molbasis[j - 3][3])) +
                                  s[j][i:i + 60])
                print("")

        # Calculate and print Mulliken Analysis (not calculated if not printed)
        if verbosity >= 2:
            print("\nMulliken Analysis")
            # First, calculate Mulliken net AO and overlap  populations
            mullikenNetAOandOvlPop = np.zeros((len(MOVectors), len(MOVectors)))
            for i in range(0, len(MOVectors)):
                for j in range(0, len(MOVectors)):
                    for k in range(0, len(MOVectors)):
                        # if, elif, else to establish the occupation number for the MO in question.
                        occ = 0
                        if ((k + 1) * 2) <= valence_electrons:
                            occ = 2
                        elif ((k + 1) * 2) > valence_electrons and (k * 2) < valence_electrons:
                            occ = 1
                        else:
                            occ = 0
                        # print(i, j, occ * (MOVectors[i][j] ** 2))
                        if i != j:
                            mullikenNetAOandOvlPop[i][j] += 2 * occ * MOVectors[i][k] * MOVectors[j][k] * overlap[i][j]
                        else:
                            mullikenNetAOandOvlPop[i][j] += occ * MOVectors[i][k] * MOVectors[j][k] * overlap[i][j]
            # Print routine for the net populations
            print("\nMuliken net AO (diagonal) and overlap (off-diagonal) populations")
            s = [""] * (len(mullikenNetAOandOvlPop) + 1)
            for i in range(0, len(mullikenNetAOandOvlPop)):
                s[0] += "{: >3}({: >3}){:>2}{}{:<2}".format(self.atoms[molbasis[i][0]].symbol,
                                                            molbasis[i][0], molbasis[i][1],
                                                            qn2symb(molbasis[i][2]),
                                                            qn2symb(molbasis[i][2], molbasis[i][3]))
                for j in range(0, len(mullikenNetAOandOvlPop)):
                    s[i + 1] += "   {: .6f} ".format(mullikenNetAOandOvlPop[i][j])
            for i in range(0, len(s[0]), 65):
                for j in range(0, (len(mullikenNetAOandOvlPop) + 1)):
                    if len(s[j]) < (i + 65):
                        if j == 0:
                            print("              " + s[j][i:len(s[0])])
                        else:
                            print("{: >3}({: >3}){:>2}{}{:<2}".format(self.atoms[molbasis[j - 1][0]].symbol,
                                                                      molbasis[j - 1][0], molbasis[j - 1][1],
                                                                      qn2symb(molbasis[j - 1][2]),
                                                                      qn2symb(molbasis[j - 1][2], molbasis[j - 1][3])) +
                                  s[j][i:len(s[0])])
                    else:
                        if j == 0:
                            print("              " + s[j][i:i + 65])
                        else:
                            print("{: >3}({: >3}){:>2}{}{:<2}".format(self.atoms[molbasis[j - 1][0]].symbol,
                                                                      molbasis[j - 1][0], molbasis[j - 1][1],
                                                                      qn2symb(molbasis[j - 1][2]),
                                                                      qn2symb(molbasis[j - 1][2], molbasis[j - 1][3])) +
                                  s[j][i:i + 65])
                print("")
            # Then, reduce Mulliken net AO and overlap populations to atoms
            mullikenRedPop = np.zeros((self.numatoms(), self.numatoms()))
            for i in range(0, len(MOVectors)):
                for j in range(0, len(MOVectors)):
                    atom1 = molbasis[i][0]
                    atom2 = molbasis[j][0]
                    mullikenRedPop[atom1][atom2] += mullikenNetAOandOvlPop[i][j]
            # Print routine for the gross populations
            print("\nMuliken reduced net atomic (diagonal) and bond (off-diagonal) populations")
            s = [""] * (len(mullikenRedPop) + 1)
            for i in range(0, len(mullikenRedPop)):
                s[0] += "   {: >3}({: >3})".format(self.atoms[i].symbol, i)
                for j in range(0, len(mullikenRedPop)):
                    s[i + 1] += " {: .6f} ".format(mullikenRedPop[i][j])
            for i in range(0, len(s[0]), 66):
                for j in range(0, (len(mullikenRedPop) + 1)):
                    if len(s[j]) < (i + 66):
                        if j == 0:
                            print("      " + s[j][i:len(s[0])])
                        else:
                            print("{: >3}({: >3})".format(self.atoms[j - 1].symbol, j - 1) + s[j][i:len(s[0])])
                    else:
                        if j == 0:
                            print("      " + s[j][i:i + 66])
                        else:
                            print("{: >3}({: >3})".format(self.atoms[j - 1].symbol, j - 1) + s[j][i:i + 66])
                print("")
            # Next calculate gross Mulliken AO populations
            mullikenGrossAOPop = np.zeros(len(MOVectors))
            for i in range(0, len(mullikenGrossAOPop)):
                for j in range(0, len(MOVectors)):
                    if i == j:
                        mullikenGrossAOPop[i] += mullikenNetAOandOvlPop[i][j]
                    else:
                        mullikenGrossAOPop[i] += mullikenNetAOandOvlPop[i][j] / 2.0
            # Print routine for the gross populations
            print("\nGross Mulliken AO populations")
            for i in range(0,len(mullikenGrossAOPop)):
                print("{: >3}({: >3}){:>2}{}{:<2} {: .6f}".format(self.atoms[molbasis[i][0]].symbol,
                                                                      molbasis[i][0], molbasis[i][1],
                                                                      qn2symb(molbasis[i][2]),
                                                                      qn2symb(molbasis[i][2], molbasis[i][3]), mullikenGrossAOPop[i]))
            print("")
            # Next calculate gross Mulliken atom populations
            mullikenGrossAtomPop = np.zeros(self.numatoms())
            for i in range(0, len(mullikenGrossAOPop)):
                mullikenGrossAtomPop[molbasis[i][0]] += mullikenGrossAOPop[i]
            print("\nGross Mulliken atomic populations")
            for i in range(0,len(mullikenGrossAtomPop)):
                print("{: >3}({: >3}) {: .6f}".format(self.atoms[i].symbol,
                                                                      i, mullikenGrossAtomPop[i]))

            print("")
            # Next determine Mulliken net atomic charges
            mullikenNetAtomCharge = np.zeros(self.numatoms())
            for i in range(0, len(mullikenNetAtomCharge)):
                mullikenNetAtomCharge[i] += self.atoms[i].valele-mullikenGrossAtomPop[i]
            print("\nNet Mulliken atomic charges")
            for i in range(0,len(mullikenGrossAtomPop)):
                print("{: >3}({: >3}) {: .6f}".format(self.atoms[i].symbol,
                                                                      i, mullikenNetAtomCharge[i]))

            print("")


        # Return the previously calculated total EHT energy
        return energy


#############################################################################################################
# Most important function so far: Read Quantum Chemistry output file and construct WellFaRe Molecule from it
#############################################################################################################

def extractCoordinates(filename, molecule, verbosity=0, distfactor=1.3, bondcutoff=0.45):
    if verbosity >= 1:
        print("\nSetting up WellFARe molecule: ", molecule.name)
    f = open(filename, 'r')
    program = "N/A"
    # Determine which QM program we're dealing with
    for line in f:
        if line.find("Entering Gaussian System, Link 0=g09") != -1:
            if verbosity >= 1:
                print("Reading Gaussian output file: ", filename)
            program = "g09"
            break
        elif line.find("* O   R   C   A *") != -1:
            if verbosity >= 1:
                print("Reading Orca output file: ", filename)
            program = "orca"
            break
    f.close()

    # GEOMETRY READING SECTION
    geom = []
    # Read through Gaussian file, read *last* "Input orientation"
    if program == "g09":
        f = open(filename, 'r')
        for line in f:
            if line.find("Input orientation:") != -1:
                if verbosity >= 2:
                    print("\nInput orientation found, reading coordinates")
                del geom[:]
                for i in range(0, 4):
                    readBuffer = f.__next__()
                while True:
                    readBuffer = f.__next__()
                    if readBuffer.find("-----------") == -1:
                        geom.append(readBuffer)
                        if verbosity >= 3:
                            readBuffer = readBuffer.split()
                            print(" Found atom: {:<3} {: .8f} {: .8f} {: .8f} in current Input orientation".format(
                                NumberToSymbol[int(readBuffer[1])], float(readBuffer[3]), float(readBuffer[4]),
                                float(readBuffer[5])))
                    else:
                        break
        if verbosity >= 1:
            print("\nReading of geometry finished.\nAdding atoms to WellFARe molecule: ", molecule.name)
        for i in geom:
            readBuffer = i.split()
            molecule.addAtom(Atom(NumberToSymbol[int(readBuffer[1])], float(readBuffer[3]), float(readBuffer[4]),
                                  float(readBuffer[5])), verbosity=verbosity)
        f.close()
    # Read through ORCA file, read *last* set of cartesian coordinates
    elif program == "orca":
        f = open(filename, 'r')
        for line in f:
            if line.find("CARTESIAN COORDINATES (ANGSTROEM)") != -1:
                if verbosity >= 2:
                    print("\nCartesian Coordinates found")
                del geom[:]
                readBuffer = f.__next__()
                while True:
                    readBuffer = f.__next__()
                    if readBuffer and readBuffer.strip():
                        geom.append(readBuffer)
                        if verbosity >= 3:
                            readBuffer = readBuffer.split()
                            print(" Found atom: {:<3} {: .8f} {: .8f} {: .8f} in current Cartesian Coordinates".format(
                                readBuffer[0], float(readBuffer[1]), float(readBuffer[2]), float(readBuffer[3])))
                    else:
                        break
        if verbosity >= 1:
            print("\nReading of geometry finished.\nAdding atoms to WellFARe molecule: ", molecule.name)
        for i in geom:
            readBuffer = i.split()
            molecule.addAtom(Atom(readBuffer[0], float(readBuffer[1]), float(readBuffer[2]), float(readBuffer[3])),
                             verbosity=verbosity)
        f.close()


        # End of routine


################################################################################
#                                                                              #
# This is the part of the program where the command line arguments are defined #
#                                                                              #
################################################################################

parser = argparse.ArgumentParser(
    description="WellFAReHMP: Wellington Fast Assessment of Reactions - Hückel MO Theory",
    epilog="recognised filetypes: g09, orca")
parser.add_argument("file", metavar='file', help="input file with structural data",
                    default="molecule.log")
# parser.add_argument("-t", "--theory", help="type of function to fit spectrum", choices=["eht", "dftb"],
#                     default="eht")
parser.add_argument("-v", "--verbosity", help="increase output verbosity", type=int, choices=[0, 1, 2, 3], default=1)

args = parser.parse_args()

###############################################################################
#                                                                             #
# The main part of the program starts here                                    #
#                                                                             #
###############################################################################

# Print GPL v3 statement and program header
ProgramHeader()

hmo_mol = Molecule("HMO Molecule")
extractCoordinates(args.file, hmo_mol, verbosity=args.verbosity)
hmo_mol.orient()

print("Total HMO Energy: {:> 16.8f}".format(hmo_mol.HMOEnergy([0], verbosity=args.verbosity)))

ProgramFooter()
