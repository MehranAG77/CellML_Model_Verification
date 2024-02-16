from sympy import symbols, Eq, simplify, Function, Add

from libcellml import Parser, Printer, Validator, cellmlElementTypeAsString

# Importing internal packages
from utilities import print_model
import cellml


variables = []      # This list contains the variables of cellml file [In CellML file

coefficients = []       # The difference between variables and coefficients is identified by the first part of their id which is like "va_12345_r1", first part indicates that this parameter is a variable

rates = []      # This list contains all the rate variables for all equations that are in the model

rate_constants = []     # This list contains the rate constants that are in the model

def equation_builder( component ):
    
    number_of_variables = component.variableCount()

    for v in range( 0, number_of_variables ):

        id = component.variable(v).id()

        identifier = id.split('_')[0]

        if identifier == 'va': variables.append( component.variable(v) )    # Since we have two different types of parameters in CellML, I put them in different lists
        elif identifier == 'co': coefficients.append( component.variable(v) )
        elif identifier == 'rc': rate_constants.append( component.variable(v) )
        elif identifier == 'ra': rates.append( component.variable(v) )

    equations=[]

    for reaction_rate in rates:
        
        rate = symbols(reaction_rate.name())

        reaction_no = reaction_rate.id().split('_')[1]
        #print(reaction_no)
        #print(len(rate_constants))
        for rate_constant in rate_constants:
            
            reaction_reactants = {}

            reaction_products = {}

            sym_list = {}

            id = rate_constant.id()
            #print(id)
            

            if id.split('_')[2] == reaction_no:

                if rate_constant.id().split('_')[1] == 'f':
                    
                    forward_rate = symbols(rate_constant.name())
                    
                elif rate_constant.id().split('_')[1] == 'r':

                    reverse_rate = symbols(rate_constant.name())

        for c_item in coefficients:

            id = c_item.id()
            #print(id)
            #print(id.split('_')[2], reaction_no)
            if id.split('_')[2] == reaction_no:

                ChEBI = id.split('_')[1]

                coefficient = c_item.name()

                value = c_item.initialValue()

                for v_item in variables:
                    #print(v_item.name())
                    if v_item.id().split('_')[1] == ChEBI:
                        variable = v_item.name()
                        #print(variable)
                        break

                if int(value) < 0:

                    reaction_reactants[variable] = coefficient

                elif int(value) > 0:

                    reaction_products[variable] = coefficient

                elif int(value) == 0:
                    print("The stoichiometric coefficient for {c} is set to zero which is wrong".format( c = coefficient ))
                    exit()

                sym_list[variable] = symbols(variable)
                sym_list[coefficient] = symbols(coefficient)

        rhs_f = forward_rate

        for item in reaction_reactants:

            rhs_f = rhs_f * sym_list[item] ** sym_list[ reaction_reactants[item] ]

        #print(rhs_f)
        rhs_r = reverse_rate

        for item in reaction_products:

            rhs_r = rhs_r * sym_list[item] ** sym_list[ reaction_products[item] ]

        #print(rhs_r)
        equations.append(Eq(rate,rhs_f-rhs_r))

    print('\nRate equations for the reaction are as below:\n                \u2193\u2193\u2193\u2193\u2193\u2193\u2193\u2193\u2193\u2193\u2193\u2193')
    for equation in equations:

            print(equation.lhs, '=', equation.rhs)


    rate_equations = []

    for v_item in variables:

        r_symbols = {}

        id = v_item.id()

        ChEBI = id.split('_')[1]

        variable = v_item.name()

        r_symbols[variable] = symbols( variable )

        constituents = {}

        for c_item in coefficients:

            c_id = c_item.id()

            c_chebi = c_id.split('_')[1]

            if c_chebi == ChEBI:

                reaction_no = c_id.split('_')[2]

                reaction_rate_coefficient = c_item.name()

                r_symbols[ reaction_rate_coefficient ] = symbols( reaction_rate_coefficient )

                for r_item in rates:

                    if r_item.id().split('_')[1] == reaction_no:

                        reaction_rate = r_item.name()

                        r_symbols[ reaction_rate ] = symbols( reaction_rate )

                        break

                constituents[reaction_rate] = reaction_rate_coefficient


        t = symbols('t')

        x = Function(r_symbols[variable])(t)

        lhs = x.diff(t)

        rhs = Add()

        for item in constituents:

            rhs += r_symbols[item] * r_symbols[constituents[item]]

        rate_equations.append(Eq(lhs,rhs))

    print('\nRate equations for variables\' concentration are as below: \n                  \u2193\u2193\u2193\u2193\u2193\u2193\u2193\u2193\u2193\u2193\u2193\u2193')
    for equation in rate_equations:

        print(equation.lhs, '=', equation.rhs)