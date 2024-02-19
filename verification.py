
"""
This function receives stoichiometric matrix (numpy array), elemental matrix (numpy array), and symbolic parameters as a rates array (sympy array)
and does some operations to check the model's validity
"""



# importing external or built-in packages
import numpy as np
import sympy as sp

def verification(stoichiometric_array, elemental_array, rate_array = 0):

    print( "\nElemental matrix is:\n", elemental_array )
    print( "\nStoichiometric matrix is:\n", stoichiometric_array )

    elemental_matrix = sp.Matrix(elemental_array)
    stoichiometric_matrix = sp.Matrix(stoichiometric_array)

    if rate_array != 0:

        conservation_matrix = elemental_matrix * stoichiometric_matrix 

        conservation_array = np.array( conservation_matrix )

        if np.all( ( conservation_array == 0 ) ) == True:

            stoichiometric_matrix_transposed =stoichiometric_matrix.transpose()
            nullspace_transposed = stoichiometric_matrix_transposed.nullspace()
            l = len(nullspace_transposed)
        
            if l == 0:
                print('There is no Left Null Space for this matrix')
            elif l == 1:
                nullspace = np.transpose(np.array(nullspace_transposed[0]))
                print('\nThe Left Null Sapce is:\n', nullspace )
                conservation_equations_array = nullspace * rate_array
                print('\nConservation equations are:\n', conservation_equations_array[0], ' = 0\n', conservation_equations_array[1], ' = 0\n' )
                return nullspace
            else:
                n_t = np.array(nullspace_transposed[0])
                counter = 1
                while counter < l:
                    n_t = np.concatenate( ( n_t, nullspace_transposed[counter] ), axis=1 )
                    counter+=1
                nullspace = np.transpose(n_t)
                print('\nThe Left Null Space is:\n', nullspace)
                conservation_equations_array = nullspace * rate_array
                print('\nConservation equations are:\n')
                
                count = 0
                while count < l:
                    print( conservation_equations_array[count], ' = 0\n')
                    count += 1
                return nullspace

        else:
            print('Conservation of Mass is violated\n')
            
            length = conservation_array.shape[1]
            for i in  range(0,length):
                if np.all(conservation_array[:, i] == 0):
                    pass
                else:
                    print ('Reaction number %i is not correctly defined' %(i+1))
                    break
            

    else:

        conservation_matrix = elemental_matrix * stoichiometric_matrix 

        conservation_array = np.array( conservation_matrix )

        if np.all( ( conservation_array == 0 ) ) == True:
        
            stoichiometric_matrix_transposed=stoichiometric_matrix.transpose()
            nullspace_transposed= stoichiometric_matrix_transposed.nullspace()
            l = len(nullspace_transposed)

            if l == 0:
                print('There is no Left Null Space for this matrix')
            elif l == 1:
                nullspace = np.transpose(np.array(nullspace_transposed[0]))
                print('\nThe Left Null Sapce is:\n', nullspace)
                return nullspace
            else:
                n_t = np.array(nullspace_transposed[0])
                counter = 1
                while counter < l:
                    n_t = np.concatenate((n_t, nullspace_transposed[counter]), axis=1)
                    counter+=1
                nullspace = np.transpose(n_t)
                print('\nThe Left Null Space is:\n', nullspace)
                conservation_equations_array = nullspace * rate_array

        else:
            print('Conservation of Mass is violated\n')
            
            length = conservation_array.shape[1]
            for i in  range(0,length):
                if np.all(conservation_array[:, i] == 0):
                    print ('Reaction number %i is not correctly defined' %i)
                    break