"""
This package gets the stoichiometric matrix, a list showing the rows of the stoichiometric matrix, and a list showing the columns of the stoichiometric matrix.
"""

import numpy as np
from sympy import symbols
import sympy as sp

import compound_element_sorter as ces

def matrix_equation_builder ( stoichiometric_matrix, rows, columns, reaction_rate_equations_dict, components ):

    """
    This function creates the concentration rate equations from the stoichiometric matrix that is constructed from the variables in CellML file
    Inputs are: rows and columns which are two dictionaries mapping row and column names to their ordering numbers
    """

    _, _ , reaction_rates, _ = ces.variable_sorter( components )

    concentration_rate_equations = {}

    for compound, row_number in rows.items():   # Since each row shows the reaction a compound participates, we go through each row and find the reactions the compound participates
        # 'rows' is a dictionary mapping compound names to their row position in the stoichiometric matrix, so we get the row number by checking the 'rows' dictionary

        temporary_reactions = {}    # To construct the right hand side of the equations, I need to store reaction name with their stoichiometric coefficient which shows the rate of consumption or production of a variable in a reaction
        # I will multiply this reaction rate to its coefficient later, so I need to keep both of them for later use as a mapping

        for column_number, element in enumerate(stoichiometric_matrix[row_number]): # now we want to get the value of the element in the stoichiometric matrix for this specific compound to see if it participates in  which reaction

            if element != 0:    # If the element value is not zero, then it participates in this reaction. We have to get the reaction name here.

                reaction_name = next((key for key, value in columns.items() if value == column_number), None)
                temporary_reactions[reaction_name] = element

        rhs = 0 # Here, I construct an empty right hand side for the equation
        
        lhs = 'd[' + compound + ']/dt'

        for reaction in temporary_reactions.keys(): # I nned to look for the reaction in the reactions list to find its name since I only have ids of the reaction which might be different with its variable name in CellML

            for reaction_rate in reaction_rates:

                if reaction == reaction_rate.id().split('_')[1]: # Here I find the CellML reaction component and retrieve its variable name
                    
                    rate_symbol = symbols(reaction_rate.name())
                    
                    rhs = rhs + temporary_reactions[reaction] * rate_symbol    # Here I need to multiply the stoichiometric coefficient element with the reaction name to construct its rate consumption equation

                    rhs = rhs.subs( rate_symbol, reaction_rate_equations_dict[reaction_rate.name()] )

        # print('\n', lhs, '=', rhs)

        concentration_rate_equations[compound] = rhs

    return concentration_rate_equations
        











            










if __name__ == '__main__':

    import CellML_reader as cmlr
    import compound_element_sorter as ces
    import equation_builder as eb

    cellml_file_dir = './docs/reactions_set.cellml'
    cellml_file = './docs/reactions_set.cellml'
    cellml_strict_mode = False

    components = cmlr.CellML_reader( cellml_file, cellml_file_dir, cellml_strict_mode )

    element_indices, compound_indices, symbols_list, compounds, stoichiometric_matrix, reaction_indices = ces.cellml_compound_element_sorter ( components )

    equations_dict = eb.equation_builder( components )

    matrix_equation_builder ( stoichiometric_matrix, compound_indices, reaction_indices, equations_dict, components )