
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





cellml_file_dir = './docs/huang_ferrell_1996.cellml'
cellml_file = './docs/huang_ferrell_1996.cellml'
cellml_strict_mode = False

components = cmlr.CellML_reader( cellml_file, cellml_file_dir, cellml_strict_mode )

reaction_rate_equations_dict= eb.equation_builder( components, print='off' )

element_indices, compound_indices, symbols_list, compounds, stoichiometric_matrix, reaction_indices = ces.cellml_compound_element_sorter ( components )

element_matrix = emb.elemental_matrix_builder( compound_indices, element_indices, compounds )

concentration_rate_equations = meb.matrix_equation_builder ( stoichiometric_matrix, compound_indices, reaction_indices, reaction_rate_equations_dict, components, print='off' )

rate_matrix = rmb.rate_matrix_builder ( symbols_list )

# Calling the function
vf.verification( stoichiometric_matrix, element_matrix, rate_matrix )

solution, time, x, sympy_to_CellML = sos.sympy_ode_solver( components, concentration_rate_equations, 40, 0.001 )

variables_to_plot = []

sos.plotter(  solution, time, variables_to_plot, x, sympy_to_CellML )