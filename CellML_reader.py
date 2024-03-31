
"""
This function receives three inputs: 1- cellml_file, which is the name of the cellml file as a string, 2- cellml_file_dir, which is the path as a string to the cellml file,
and 3- cellml_strict_mode, which is a boolean with True and False identifying CellML 2.0 and 1.1 versions, respectively.
The new version (updated on 12/12/2023) returns the class of a component (libcellml class) that can be used to find variables and coefficients by another module
"""


# Importing external packages
from libcellml import Parser, Printer, Validator, cellmlElementTypeAsString
import numpy as np
import chemparse as chp

# Importing internal packages
from utilities import print_model
import cellml
import chebi_fetch as chf
import sys



cellml_file_dir = './docs/huang_ferrell_1996.cellml'
cellml_file = './docs/huang_ferrell_1996.cellml'
cellml_strict_mode = False

def CellML_reader( cellml_file, cellml_file_dir, cellml_strict_mode ):        # strict_mode is True when CellML version is 2.0 and False when 1.0 and 1.1

    model = cellml.parse_model(cellml_file, cellml_strict_mode)
    if cellml.validate_model(model) > 0:
        exit(-1)

    importer = cellml.resolve_imports(model, cellml_file_dir, cellml_strict_mode)
    if model.hasUnresolvedImports():
        print("unresolved imports?")
        exit(-2)

    if cellml.validate_model(model) > 0:
        print('Validation issues found')
        exit(-3)

    print('Model was parsed, resolved, and validated without any issues.')

    # need a flattened model for analysing
    flat_model = cellml.flatten_model(model, importer)
    if cellml.validate_model(flat_model) > 0:
        print('Validation issues found in flattened model')
        exit(-4)

    print('Model was flattened without any issues.')

    

    number_of_components = model.componentCount()
    #print('numer of components:  ', number_of_components)
    
    if number_of_components >= 1 :

        components = []

        for component in range(0, number_of_components):

            components.append( model.component( component ) )

    else:
        print( 'There is no component in the CellML file' )
        exit(-4)

    return components


if __name__ == '__main__':

    print(CellML_reader( cellml_file, cellml_file_dir, cellml_strict_mode ))