"""This module contains a range of tools for evaluating protein structure."""

from .contact_order import calculate_contact_order
from .hydrophobic_fitness import calculate_hydrophobic_fitness
from .dssp import tag_dssp_data
from .pacc import fit_heptad_register, PACCAnalysis

del contact_order
del hydrophobic_fitness
