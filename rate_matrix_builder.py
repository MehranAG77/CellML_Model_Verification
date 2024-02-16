"""
This function receives a list containing the compounds as strings and returns a sympy matrix containing the parameters of the compounds
"""


# Importing external or built-in packages
import numpy as np
import sympy as sp


def rate_mat_builder( coms ):

    symbols_list = sp.symbols(coms)

    rate_mat = sp.Matrix(symbols_list)

    return rate_mat
        