
"""
This function receives stoichiometric matrix (numpy array), elemental matrix (numpy array), and symbolic parameters as a rates array (sympy array)
and does some operations to check the model's validity
"""



# importing external or built-in packages
import numpy as np
import sympy as sp

def verification(st_array, el_array, rate_array = 0):

    print( "\nElemental matrix is:\n", el_array )
    print( "\nStoichiometric matrix is:\n", st_array )

    el_array = sp.Matrix(el_array)
    st_array = sp.Matrix(st_array)

    if rate_array != 0:

        con_mat = el_array * st_array

        con_mat = np.array( con_mat )

        if np.all( ( con_mat == 0 ) ) == True:

            st_array_tr =st_array.transpose()
            tr_null = st_array_tr.nullspace()
            l = len(tr_null)
        
            if l == 0:
                print('There is no Left Null Space for this matrix')
            elif l == 1:
                n = np.transpose(np.array(tr_null[0]))
                print('\nThe Left Null Sapce is:\n', n)
                con_eqs_arr = n * rate_array
                print('\nConservation equations are:\n', con_eqs_arr[0], ' = 0\n', con_eqs_arr[1], ' = 0\n' )
                return n
            else:
                n_t = np.array(tr_null[0])
                counter = 1
                while counter < l:
                    n_t = np.concatenate((n_t, tr_null[counter]), axis=1)
                    counter+=1
                n = np.transpose(n_t)
                print('\nThe Left Null Space is:\n', n)
                con_eqs_arr = n * rate_array
                print('\nConservation equations are:\n')
                
                count = 0
                while count < l:
                    print( con_eqs_arr[count], ' = 0\n')
                    count += 1
                return n

        else:
            print('Conservation of Mass is violated\n')
            
            length = con_mat.shape[1]
            for i in  range(0,length):
                if np.all(con_mat[:, i] == 0):
                    pass
                else:
                    print ('Reaction number %i is not correctly defined' %(i+1))
                    break
            

    else:

        con_mat = el_array * st_array

        con_mat = np.array( con_mat )

        if np.all( ( con_mat == 0 ) ) == True:
        
            st_array_tr=st_array.transpose()
            tr_null= st_array_tr.nullspace()
            l = len(tr_null)

            if l == 0:
                print('There is no Left Null Space for this matrix')
            elif l == 1:
                n = np.transpose(np.array(tr_null[0]))
                print('\nThe Left Null Sapce is:\n', n)
                return n
            else:
                n_t = np.array(tr_null[0])
                counter = 1
                while counter < l:
                    n_t = np.concatenate((n_t, tr_null[counter]), axis=1)
                    counter+=1
                n = np.transpose(n_t)
                print('\nThe Left Null Space is:\n', n)
                con_eqs_arr = n * rate_array

        else:
            print('Conservation of Mass is violated\n')
            
            length = con_mat.shape[1]
            for i in  range(0,length):
                if np.all(con_mat[:, i] == 0):
                    print ('Reaction number %i is not correctly defined' %i)
                    break