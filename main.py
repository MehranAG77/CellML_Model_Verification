
"""
    *******************************************************************************************************************************
    *                                                                                                                             *
    *  This module reads CellML files containing reaction stoichiometries and equations                                           *
    *                                                                                                                             *
    *  Then using the modules provided constructs stoichiometric and elemental matrices and verifies the model                    *
    *                                                                                                                             *
    *  Then the corresponding Concentration Rate Equations are generated using the Stoichiometric matrix and solved               *
    *                                                                                                                             *
    *******************************************************************************************************************************
"""

# Importing external or built-in packages
import numpy as np
import sympy as sp
import os

# Importing internal packages
import excel_read as xls
import compounds_extractor as comex
import stochiometric_matrix_builer as smb
import elemental_matrix_builder as emb
import rate_matrix_builder as rmb
import verification as vf
import CellML_reader as cmlr
import compound_element_sorter as ces
import equation_builder as eb
import matrix_equation_builder as meb
import sympy_ode_solver as sos

command = 'cls' if os.name == 'nt' else 'clear'
os.system(command)



cellml_file_dir = './docs/reactions_set.cellml'
cellml_file = './docs/reactions_set.cellml'
cellml_strict_mode = False

components = cmlr.CellML_reader( cellml_file, cellml_file_dir, cellml_strict_mode )

variables, coefficients, rates, rate_constants, boundary_conditions = ces.variable_sorter( components )

reaction_rate_equations_dict, bc_equations = eb.equation_builder( components, 'on' ) #print

element_indices, compound_indices, reaction_indices, symbols_list, compounds, bc_coefficients = ces.cellml_compound_element_sorter ( components )

element_matrix = emb.elemental_matrix_builder( compound_indices, element_indices, compounds )

stoichiometric_matrix = smb.stoichiometric_matrix_builder( reaction_indices, compound_indices, coefficients, bc_coefficients )

concentration_rate_equations = meb.matrix_equation_builder ( stoichiometric_matrix, compound_indices, reaction_indices, reaction_rate_equations_dict, components, 'on' ) #print

rate_matrix = rmb.rate_matrix_builder ( symbols_list )

# Calling the function
vf.verification( stoichiometric_matrix, element_matrix, rate_matrix )

solution, time, x, sympy_to_CellML = sos.sympy_ode_solver( components, concentration_rate_equations, 40, 0.001 )

variables_to_plot = []

sos.plotter(  solution, time, variables_to_plot, x, sympy_to_CellML )