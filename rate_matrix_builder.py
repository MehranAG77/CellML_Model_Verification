"""
This function receives a list containing the compounds as strings and returns a sympy matrix containing the parameters of the compounds
"""


# Importing external or built-in packages
import numpy as np
import sympy as sp


def rate_matrix_builder( compounds ):

    symbols_list = sp.symbols(compounds)

    rate_matrix = sp.Matrix(symbols_list)

    return rate_matrix
        