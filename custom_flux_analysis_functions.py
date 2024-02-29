import numpy as np
import numpy.linalg as linalg
from sympy import Matrix
from scipy.linalg import null_space

# print_reactions( N, v_lab, m_lab, rev )
# based on the stoichiometric matrix and cell arrays containing the labes for colmunms and
# rows print the network reactions on the screen
# Input arguments:
# N - mxn matrix containing the stoichiometric coeficients for the m metabolites and n reactions
# v_lab - cell array containing the reaction names
# m_lab - cell array containing the metabolite names
# Optional arguments:
# rev - array with logical entries (0,1) describing 1-reversible, 0-irreversible



# Functions created by Prof. Dr. Sebastian Aljoscha Wahl
def print_reactions(N, v_lab, m_lab, revv=0):
    if revv == 0:
        rev = np.zeros(N.shape[1])
    else:
        rev = revv
    
    # go through all columns and print the reactions
    print('Reactions of the network:')

    for i in range(N.shape[1]):
        # positive stoichiometric coefficient == product of the reaction
        p = np.where(N[:,i] > 0)
        p = p[0]
        # negative stoichiometric coefficient == substrate of the reaction
        s = np.where(N[:,i] < 0)
        s = s[0]

        # we are at column i, print out the name of the reaction
        print('%g %s:' % (i+1, v_lab[i]), end='')

        if len(s) != 0: # is there at least one substrate?
            # left side = substrates
            print('    %g %s ' % (-N[s[0], i], m_lab[s[0]]), end='') # print first substrate
            for ii in range(1, len(s)): # are there more?
                print('+ %g %s ' % (-N[s[ii], i], m_lab[s[ii]]), end='')
        else:
            print('   ', end='')
    
        if rev[i] > 0:
            print(' <--> ', end='')
        else:
            print(' ---> ', end='')
        
        if len(p) != 0: # is there at least one product?
            # right side = product
            if len(p) == 1: # new line needed when there is only one product
                print(' %g %s' % (N[p[0], i], m_lab[p[0]])) # print first product
            else:
                print(' %g %s ' % (N[p[0], i], m_lab[p[0]]), end='') # print first product
                for ii in range(1, len(p)): # are there more?
                    print(' + %g %s ' % (N[p[ii], i], m_lab[p[ii]]))
        else:
            print('')
          
        

# Customized 'print_reactions' to transfer the reaction equations into the desired csv form
def print_reactions_csv(N, v_lab, m_lab, flux,ex_flx):
    
    # go through all columns and print the reactions
    print('-'*12,'CSV-reaction form of the network','-'*12,'\n')
    print('Reaction Formula,Relative Flux,Reaction Name')
    
    # sort out reactions that are smaller than 0.1% of the uptake rate
    index = []
    for n,m in enumerate(flux):
        if m >= 0.001 and ex_flx[n] == 0:
            index.append(n)
    
    for i in index:
        # positive stoichiometric coefficient == product of the reaction
        p = np.where(N[:,i] > 0)
        p = p[0]
        # negative stoichiometric coefficient == substrate of the reaction
        s = np.where(N[:,i] < 0)
        s = s[0]


        if len(s) != 0: # is there at least one substrate?
            # left side = substrates
            print('%g %s ' % (-N[s[0], i], m_lab[s[0]]), end='') # print first substrate
            for ii in range(1, len(s)): # are there more?
                print('+ %g %s ' % (-N[s[ii], i], m_lab[s[ii]]), end='')
        else:
            print('   ', end='')
    
        print(' <=> ', end='')
        
        if len(p) != 0: # is there at least one product?
            # right side = product
            if len(p) == 1: # new line needed when there is only one product
                print(' %g %s ' % (N[p[0], i], m_lab[p[0]]), end='') # print first product
            else:
                print(' %g %s ' % (N[p[0], i], m_lab[p[0]]), end='') # print first product
                for ii in range(1, len(p)): # are there more?
                    print('+ %g %s ' % (N[p[ii], i], m_lab[p[ii]]),end='')
        else:
            print('', end='')
            
        # print the optimal flux and name of the reaction into the csv-file
        print( ',{:0.4f},{:s}'.format(flux[i], v_lab[i]))


# Functions created by Prof. Dr. Sebastian Aljoscha Wahl        
def check_blocked_reactions(N, v_lab):
    print('-'*12,'Blocked reactions?','-'*12)
    NS = null_space(N)
    n_blocked = []

    for i in range(NS.shape[0]):
        if all(np.absolute(NS[i,:]) < 1e-12): # entries < 1e-12?
            print('Reaction %s is blocked' % (v_lab[i]))
            n_blocked = np.append(n_blocked, i)
        
    if len(n_blocked) == 0:
        print('There is no blocked reactions in this metabolic network')
    
    return n_blocked