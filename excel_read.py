

"""
This function reads an Excel file containing two rows of headers: First row containing the name of the reaction compound
and The Second row containing an alphabet of 'E' or 'C' as a sign of 'Element' or 'Coefficient', respectively.
And returns a datafrmae containing the reaction compound elements.
So the headers =[0,1] indicates that rows 0 and 1 should be treated as headers.
"""


# Importing external or built-in packages
import numpy as np
import pandas as pd

def read_excel( path , headers = [0,1] ):

    df = pd.read_excel(path, header=headers)

    return df




"""
This function reads an Excel file having The First row as a header and the rest as variables
and returns a numpy array containing the compound elements and their respective subscript
"""

def read_excel_reactions ( path ):

    df = pd.read_excel(path)

    dfn = df.to_numpy()

    return dfn