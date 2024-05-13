from sympy import symbols, lambdify
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
from compound_element_sorter import variable_name_mapper

def sympy_ode_solver( components, concentration_rate_equations ):

    """
    This function solves a set of Ordinary Differential Equations (ODEs) by Sympy internal solver which is Runge Kutta 4th order.
    It gets a dictionary mapping the 'compound', represented by its chemical formula, to the equation created by Sympy variables using the stoichiometric matrix
    'components' is the list containing the imported CellML components
    """

    all_symbols = set().union(*[eq.free_symbols for eq in concentration_rate_equations.values()])   # I need to know all variables that I have in these equations. This is why I should have replaced all variables that have values and are not meant ot be solved for

    number_of_variables = str(len(all_symbols)) # I need to count the number of variables that I have so I can create the list of sympy variables as 'x' shown below

    sympy_to_CellML = {}    # I need to map sympy variables to CellML variables written by user. This is because I have ChEBI code and variable name for the variables that I have
                            # For construction of the stoichiometric matrix, I used ChEBI codes and here I need variable names

    sympy_equations ={}

    chebi_to_CellML, chebi_initial_values = variable_name_mapper( components )

    initial_values = [] # Initial values should be stored in the same order as variables in 'x'

    xdot = []   # Sympy needs the equations in the set as a list containing all the equations

    x = symbols( 'x:' + number_of_variables )# Creation of the number of Sympy variables that I need

    # I need to replace the CellML variables I have used to write the equations with the Sympy variables which are like 'x0', 'x1', 'x3', ...

    for i, variable in enumerate ( all_symbols ):

        # here I will map the sympy variable to the CellML variable and store them in a dictionary for future reference
        sympy_to_CellML[x[i]] = str(variable)

        # Now, I need to replace the CellML variables with the corresponding Sympy variable in all equations of concentration rates
        # So, I need a 'for loop' to check all the equations
        for key in concentration_rate_equations.keys():

            concentration_rate_equations[key] =concentration_rate_equations[key].subs( variable, x[i] )

    for key in concentration_rate_equations.keys():

        # key in the concentration_rate_equations dictionary is the compound name as a string in chemical formula format

        to_find = chebi_to_CellML[key]  # I get the CellML variable of the chemical formulation. then I need to find its correspondent Sympy variable. Hence, I need a 'for' loop as below

        for x_find, value in sympy_to_CellML.items():
            
            # to_find is a CellML variable, x_find is sympy variable, value is again CellMl variable mapped from x_find. So, I'll compare value and to_find
            #If they are the same, I'll map the sympy variable to its corresponding equation in a new dictionary which belongs to it

            if value == to_find:

                sympy_equations[x_find] = concentration_rate_equations[key]
                break

    # Now, I have the equations all in sympy format. A Sympy type variable mapped to its corresponding equation in Sympy formatted variables
    # Now, I have to find the initial values for these Sympy variables to solve the set of ODEs
    for i in range( len(x) ):

        variable = sympy_to_CellML[symbols('x'+str(i))]

        for key, value in chebi_to_CellML.items():

            if variable == value:

                initial_values.append( float(chebi_initial_values[key]) )
                break
        
        variable = symbols('x'+str(i))

        # xdot is the list of equations that Sympy needs to solve
        for key in sympy_equations.keys():

            if variable == key:

                xdot.append(sympy_equations[key])
                break


    #################################################################
    #################################################################

    # Now, lets solve the obtained equations

    # The time step and End point for the solution
    t_max = 10
    delta_t = 0.001




    t = symbols( 't' )
    f = lambdify( ( t, x ), xdot )

    n = int(t_max/delta_t)+1

    t_eval = np.linspace( 0, t_max, n )

    solution = scipy.integrate.solve_ivp( f, (0, t_max), initial_values, t_eval = t_eval )

    y = solution.y

    legends = []


    # This 'for' loop creates the legend based on the CellMl variable names for the independent variables in ODEs
    for variable in x:

        legends.append( sympy_to_CellML[variable] )

    plt.plot( t_eval, y.T )
    plt.title( 'Chemical Kinetics' )
    plt.legend( legends, shadow = True )
    plt.xlabel('time')
    plt.ylabel('concentration')

    plt.show()





if __name__ == '__main__':

    import CellML_reader as cmlr
    import compound_element_sorter as ces
    import equation_builder as eb
    import matrix_equation_builder as meb

    cellml_file_dir = './docs/NitrosylBromide.cellml'
    cellml_file = './docs/NitrosylBromide.cellml'
    cellml_strict_mode = False

    components = cmlr.CellML_reader( cellml_file, cellml_file_dir, cellml_strict_mode )

    element_indices, compound_indices, symbols_list, compounds, stoichiometric_matrix, reaction_indices = ces.cellml_compound_element_sorter ( components )

    equations_dict = eb.equation_builder( components )

    equations = meb.matrix_equation_builder ( stoichiometric_matrix, compound_indices, reaction_indices, equations_dict, components )

    sympy_ode_solver( components, equations )