#!/usr/bin/env python
# coding: utf-8


import os
import sys, getopt
import numpy as np
import xgboost as xgb

from collections import namedtuple
from enum import Enum


class AAHydrophobicityClass(Enum):
    NEUTRAL = 0
    HYDROPHILIC = 1
    HYDROPHOBIC = 2


class AAChemicalProperty(Enum):
    BASIC = 0
    AMIDE = 1
    ACIDIC = 2
    SULFUR = 3
    HYDROXYL = 4
    AROMATIC = 5
    ALIPHATIC = 6


class AASize(Enum):
    SMALL = 0
    MEDIUM = 1
    LARGE = 2
    VERY_LARGE = 3
    VERY_SMALL = 4


class AAHydrogenBonding(Enum):
    NONE = 0
    DONOR = 1
    DONOR_ACCEPTOR = 2
    ACCEPTOR = 3


class AAPolarity(Enum):
    NONPOLAR = 0
    POLAR_BASIC = 1
    POLAR_NEUTRAL = 2
    POLAR_ACIDIC = 3


AAProperty = namedtuple(
    "AAProperty",
    "name_3_letter, name_1_letter, volume, hydrophobicity_value,"
    "flexibility, chemical_property, size, "
    "polarity, hydrogen_bonding, hydrophobicity_class",
)

aa3_map = {}

aa3_map["ALA"] = AAProperty(
    name_3_letter="ALA",
    name_1_letter="A",
    volume=88.6,
    hydrophobicity_value=0,
    flexibility=1,
    chemical_property=AAChemicalProperty.ALIPHATIC,
    size=AASize.VERY_SMALL,
    polarity=AAPolarity.NONPOLAR,
    hydrogen_bonding=AAHydrogenBonding.NONE,
    hydrophobicity_class=AAHydrophobicityClass.HYDROPHOBIC,
)

aa3_map["CYS"] = AAProperty(
    name_3_letter="CYS",
    name_1_letter="C",
    volume=108.5,
    hydrophobicity_value=0.49,
    flexibility=3,
    chemical_property=AAChemicalProperty.SULFUR,
    size=AASize.SMALL,
    polarity=AAPolarity.NONPOLAR,
    hydrogen_bonding=AAHydrogenBonding.NONE,
    hydrophobicity_class=AAHydrophobicityClass.HYDROPHOBIC,
)

aa3_map["ASP"] = AAProperty(
    name_3_letter="ASP",
    name_1_letter="D",
    volume=111.1,
    hydrophobicity_value=2.95,
    flexibility=18,
    chemical_property=AAChemicalProperty.ACIDIC,
    size=AASize.SMALL,
    polarity=AAPolarity.POLAR_ACIDIC,
    hydrogen_bonding=AAHydrogenBonding.ACCEPTOR,
    hydrophobicity_class=AAHydrophobicityClass.HYDROPHILIC,
)

aa3_map["GLU"] = AAProperty(
    name_3_letter="GLU",
    name_1_letter="E",
    volume=138.4,
    hydrophobicity_value=1.64,
    flexibility=54,
    chemical_property=AAChemicalProperty.ACIDIC,
    size=AASize.MEDIUM,
    polarity=AAPolarity.POLAR_ACIDIC,
    hydrogen_bonding=AAHydrogenBonding.ACCEPTOR,
    hydrophobicity_class=AAHydrophobicityClass.HYDROPHILIC,
)

aa3_map["PHE"] = AAProperty(
    name_3_letter="PHE",
    name_1_letter="F",
    volume=189.9,
    hydrophobicity_value=-2.2,
    flexibility=18,
    chemical_property=AAChemicalProperty.AROMATIC,
    size=AASize.VERY_LARGE,
    polarity=AAPolarity.NONPOLAR,
    hydrogen_bonding=AAHydrogenBonding.NONE,
    hydrophobicity_class=AAHydrophobicityClass.HYDROPHOBIC,
)

aa3_map["GLY"] = AAProperty(
    name_3_letter="GLY",
    name_1_letter="G",
    volume=60.1,
    hydrophobicity_value=1.72,
    flexibility=1,
    chemical_property=AAChemicalProperty.ALIPHATIC,
    size=AASize.VERY_SMALL,
    polarity=AAPolarity.NONPOLAR,
    hydrogen_bonding=AAHydrogenBonding.NONE,
    hydrophobicity_class=AAHydrophobicityClass.NEUTRAL,
)

aa3_map["HIS"] = AAProperty(
    name_3_letter="HIS",
    name_1_letter="H",
    volume=153.2,
    hydrophobicity_value=4.76,
    flexibility=36,
    chemical_property=AAChemicalProperty.BASIC,
    size=AASize.MEDIUM,
    polarity=AAPolarity.POLAR_BASIC,
    hydrogen_bonding=AAHydrogenBonding.DONOR_ACCEPTOR,
    hydrophobicity_class=AAHydrophobicityClass.NEUTRAL,
)

aa3_map["ILE"] = AAProperty(
    name_3_letter="ILE",
    name_1_letter="I",
    volume=166.7,
    hydrophobicity_value=-1.56,
    flexibility=9,
    chemical_property=AAChemicalProperty.ALIPHATIC,
    size=AASize.LARGE,
    polarity=AAPolarity.NONPOLAR,
    hydrogen_bonding=AAHydrogenBonding.NONE,
    hydrophobicity_class=AAHydrophobicityClass.HYDROPHOBIC,
)

aa3_map["LYS"] = AAProperty(
    name_3_letter="LYS",
    name_1_letter="K",
    volume=168.6,
    hydrophobicity_value=5.39,
    flexibility=81,
    chemical_property=AAChemicalProperty.BASIC,
    size=AASize.LARGE,
    polarity=AAPolarity.POLAR_BASIC,
    hydrogen_bonding=AAHydrogenBonding.DONOR,
    hydrophobicity_class=AAHydrophobicityClass.HYDROPHILIC,
)

aa3_map["LEU"] = AAProperty(
    name_3_letter="LEU",
    name_1_letter="L",
    volume=166.7,
    hydrophobicity_value=-1.81,
    flexibility=9,
    chemical_property=AAChemicalProperty.ALIPHATIC,
    size=AASize.LARGE,
    polarity=AAPolarity.NONPOLAR,
    hydrogen_bonding=AAHydrogenBonding.NONE,
    hydrophobicity_class=AAHydrophobicityClass.HYDROPHOBIC,
)

aa3_map["MET"] = AAProperty(
    name_3_letter="MET",
    name_1_letter="M",
    volume=162.9,
    hydrophobicity_value=-0.76,
    flexibility=27,
    chemical_property=AAChemicalProperty.SULFUR,
    size=AASize.LARGE,
    polarity=AAPolarity.NONPOLAR,
    hydrogen_bonding=AAHydrogenBonding.NONE,
    hydrophobicity_class=AAHydrophobicityClass.HYDROPHOBIC,
)

aa3_map["ASN"] = AAProperty(
    name_3_letter="ASN",
    name_1_letter="N",
    volume=114.1,
    hydrophobicity_value=3.47,
    flexibility=36,
    chemical_property=AAChemicalProperty.AMIDE,
    size=AASize.SMALL,
    polarity=AAPolarity.POLAR_NEUTRAL,
    hydrogen_bonding=AAHydrogenBonding.DONOR_ACCEPTOR,
    hydrophobicity_class=AAHydrophobicityClass.HYDROPHILIC,
)

aa3_map["PRO"] = AAProperty(
    name_3_letter="PRO",
    name_1_letter="P",
    volume=112.7,
    hydrophobicity_value=-1.52,
    flexibility=2,
    chemical_property=AAChemicalProperty.ALIPHATIC,
    size=AASize.SMALL,
    polarity=AAPolarity.NONPOLAR,
    hydrogen_bonding=AAHydrogenBonding.NONE,
    hydrophobicity_class=AAHydrophobicityClass.NEUTRAL,
)

aa3_map["GLN"] = AAProperty(
    name_3_letter="GLN",
    name_1_letter="Q",
    volume=143.8,
    hydrophobicity_value=3.01,
    flexibility=108,
    chemical_property=AAChemicalProperty.AMIDE,
    size=AASize.MEDIUM,
    polarity=AAPolarity.POLAR_NEUTRAL,
    hydrogen_bonding=AAHydrogenBonding.DONOR_ACCEPTOR,
    hydrophobicity_class=AAHydrophobicityClass.HYDROPHILIC,
)

aa3_map["ARG"] = AAProperty(
    name_3_letter="ARG",
    name_1_letter="R",
    volume=173.4,
    hydrophobicity_value=3.71,
    flexibility=81,
    chemical_property=AAChemicalProperty.BASIC,
    size=AASize.LARGE,
    polarity=AAPolarity.POLAR_BASIC,
    hydrogen_bonding=AAHydrogenBonding.DONOR,
    hydrophobicity_class=AAHydrophobicityClass.HYDROPHILIC,
)

aa3_map["SER"] = AAProperty(
    name_3_letter="SER",
    name_1_letter="S",
    volume=89.0,
    hydrophobicity_value=1.83,
    flexibility=3,
    chemical_property=AAChemicalProperty.HYDROXYL,
    size=AASize.VERY_SMALL,
    polarity=AAPolarity.POLAR_NEUTRAL,
    hydrogen_bonding=AAHydrogenBonding.DONOR_ACCEPTOR,
    hydrophobicity_class=AAHydrophobicityClass.NEUTRAL,
)

aa3_map["THR"] = AAProperty(
    name_3_letter="THR",
    name_1_letter="T",
    volume=116.1,
    hydrophobicity_value=1.78,
    flexibility=3,
    chemical_property=AAChemicalProperty.HYDROXYL,
    size=AASize.SMALL,
    polarity=AAPolarity.POLAR_NEUTRAL,
    hydrogen_bonding=AAHydrogenBonding.DONOR_ACCEPTOR,
    hydrophobicity_class=AAHydrophobicityClass.NEUTRAL,
)

aa3_map["VAL"] = AAProperty(
    name_3_letter="VAL",
    name_1_letter="V",
    volume=140.0,
    hydrophobicity_value=-0.78,
    flexibility=3,
    chemical_property=AAChemicalProperty.ALIPHATIC,
    size=AASize.MEDIUM,
    polarity=AAPolarity.NONPOLAR,
    hydrogen_bonding=AAHydrogenBonding.NONE,
    hydrophobicity_class=AAHydrophobicityClass.HYDROPHOBIC,
)

aa3_map["TRP"] = AAProperty(
    name_3_letter="TRP",
    name_1_letter="W",
    volume=227.8,
    hydrophobicity_value=-0.38,
    flexibility=36,
    chemical_property=AAChemicalProperty.AROMATIC,
    size=AASize.VERY_LARGE,
    polarity=AAPolarity.NONPOLAR,
    hydrogen_bonding=AAHydrogenBonding.DONOR,
    hydrophobicity_class=AAHydrophobicityClass.HYDROPHOBIC,
)

aa3_map["TYR"] = AAProperty(
    name_3_letter="TYR",
    name_1_letter="Y",
    volume=193.6,
    hydrophobicity_value=-1.09,
    flexibility=18,
    chemical_property=AAChemicalProperty.AROMATIC,
    size=AASize.VERY_LARGE,
    polarity=AAPolarity.POLAR_NEUTRAL,
    hydrogen_bonding=AAHydrogenBonding.DONOR_ACCEPTOR,
    hydrophobicity_class=AAHydrophobicityClass.NEUTRAL,
)


aa1_map = {aa.name_1_letter: aa for k, aa in aa3_map.items()}

aa3to1 = {aa.name_3_letter: aa.name_1_letter for k, aa in aa1_map.items()}
aa1to3 = {aa.name_1_letter: aa.name_3_letter for k, aa in aa1_map.items()}

aa3_label = {aa: f"{i+1}" for i, aa in enumerate(sorted(aa3to1.keys()))}
aa3_label["GLU"], aa3_label["GLN"] = aa3_label["GLN"], aa3_label["GLU"]

aa1_label = {aa3_map[aa].name_1_letter: f"{i}" for aa, i in aa3_label.items()}
aa1_mttype_label = {a: str(i + 1) for i, a in enumerate("AFCDNEQGHLIKMPRSTVWY")}


def net_flexibility2(aa1_map, wild, mutant):
    wt_flex = aa1_map[wild].flexibility if wild in aa1_map else 0
    mt_flex = aa1_map[mutant].flexibility if mutant in aa1_map else 0
    return int(mt_flex) - int(wt_flex)


def net_volume2(aa1_map, wild, mutant):
    wt_volm = aa1_map[wild].volume if wild in aa1_map else 0
    mt_volm = aa1_map[mutant].volume if mutant in aa1_map else 0
    return "{:0.1f}".format(float(mt_volm) - float(wt_volm))


def net_hydrophobicity2(aa1_map, wild, mutant):
    """Return ddG(mt-wt)

    Parameters
    ----------
    wild: str 1 letter code of the wild-type amino acid
    mutation: str 1 letter code of the mutant amino acid
    Reference: https://doi.org/10.1073/pnas.1103979108
    """
    wt_hydphob = aa1_map[wild].hydrophobicity_value if wild in aa1_map else 0
    mt_hydphob = aa1_map[mutant].hydrophobicity_value if mutant in aa1_map else 0
    return "{:0.1f}".format(float(mt_hydphob) - float(wt_hydphob))


def mutation_chemical2(aa1_map, wild, mutant):
    wt_chem = aa1_map[wild].chemical_property.value
    mt_chem = aa1_map[mutant].chemical_property.value
    n_chem_props = len(AAChemicalProperty)
    label_chem = (wt_chem + 1) % n_chem_props * n_chem_props + mt_chem
    return label_chem


def mutation_hydrophobicity2(aa1_map, wild, mutant):
    wt_hydrphbcls = aa1_map[wild].hydrophobicity_class.value
    mt_hydrphbcls = aa1_map[mutant].hydrophobicity_class.value
    n_hydrphbcls = len(AAHydrophobicityClass)
    label_indx = (wt_hydrphbcls + 1) % n_hydrphbcls * n_hydrphbcls + mt_hydrphbcls
    return label_indx


def mutation_polarity2(aa1_map, wild, mutant):
    wt_pol = aa1_map[wild].polarity.value
    # Done to emulate the manual label-encoding used in original code
    if aa1_map[wild].polarity in [AAPolarity.NONPOLAR, AAPolarity.POLAR_BASIC]:
        wt_pol = 1 - wt_pol
    mt_pol = aa1_map[mutant].polarity.value
    n_polcls = len(AAPolarity)
    label_polarity = (wt_pol) % n_polcls * n_polcls + mt_pol
    return label_polarity


def mutation_size2(aa1_map, wild, mutant):
    wt_size = aa1_map[wild].size.value
    mt_size = aa1_map[mutant].size.value
    n_size = len(AASize)
    label_size = (wt_size + 1) % n_size * n_size + mt_size
    return label_size


def mutation_hbonds2(aa1_map, wild, mutant):
    wt_hb = aa1_map[wild].hydrogen_bonding.value
    # Done to emulate the manual label-encoding used in original code
    if aa1_map[wild].hydrogen_bonding in [
        AAHydrogenBonding.NONE,
        AAHydrogenBonding.DONOR,
    ]:
        wt_hb = 1 - wt_hb
    mt_hb = aa1_map[mutant].hydrogen_bonding.value
    n_hb = len(AAHydrogenBonding)
    label_hb = (wt_hb) % n_hb * n_hb + mt_hb
    return label_hb


def mutation_type2(aa1_label, wild, mutant):
    wt_lbl, mt_lbl = int(aa1_label.get(wild, 0)), int(aa1_label.get(mutant, 0))
    if wild == mutant:
        return (wt_lbl, mt_lbl, 0)
    if mt_lbl > wt_lbl:
        mt_lbl = mt_lbl - 1
    mutation_lbl = (int(wt_lbl) - 1) * 19 + (int(mt_lbl))
    return mutation_lbl
