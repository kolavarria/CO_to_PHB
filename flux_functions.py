#copied from M2_functions by Aljoscha Wahl
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

def print_reactions(N, v_lab, m_lab,flux):#revv=0
    #if revv == 0:
    #    rev = np.zeros(N.shape[1])
    #else:
     #   rev = revv
    
    # go through all columns and print the reactions
    print('Reactions of the network:')
    
    for i in range(N.shape[1]):
        # positive stoichiometric coefficient == substrate of the reaction
        p = np.where(N[:,i] > 0)
        p = p[0]
        # negative stoichiometric coefficient == product of the reaction
        s = np.where(N[:,i] < 0)
        s = s[0]
        
        if flux[i]<= 3e-9:
            flux[i] = 0

        # we are at column i, print out the name of the reaction
        #print('%s:'% (v_lab[i]), end='')

        if len(s) != 0: # is there at least one substrate?
            # left side = substrates
            if len(s) == 1:
                if -N[s[0],i]==1:
                    print('%s'% (m_lab[s[0]]), end='') # print first substrate
                else:
                    print('%g %s' % (-N[s[0], i], m_lab[s[0]]), end='')
            else:
                if -N[s[0],i]==1:
                    print('%s'% (m_lab[s[0]]),end='')
                else:
                    print('%g %s'%(-N[s[0], i], m_lab[s[0]]), end='')
                for ii in range(1, len(s)): # are there more?
                    if -N[s[ii],i]==1:
                        print(' + %s'%(m_lab[s[ii]]), end='')
                    else:
                        print(' + %g %s'% (-N[s[ii], i], m_lab[s[ii]]), end='')
                #print('<=>',end='')
        else:
            print('', end='')
    
        #if rev[i] > 0:
         #   print('<=>', end='')
        #else:
         #   print('->', end='')
        print(' <=> ', end='')   
        
        if len(p) != 0: # is there at least one product?
            # right side = product
            if len(p) == 1: # new line needed when there is only one product
                if N[p[0],i]==1:
                    print('%s'% (m_lab[p[0]]), end='')# print first product
            else:
                if N[p[0],i]==1:
                    print('%s'% ( m_lab[p[0]]),end='')# print first product
                else:
                    print('%g %s'%(N[p[0], i], m_lab[p[0]]),end='')
                for ii in range(1, len(p)): # are there more?
                        if N[p[ii],i]==1: 
                            print(' + %s'% (m_lab[p[ii]]),end='')
                        else:
                            print(' + %g %s'%(N[p[ii], i], m_lab[p[ii]]),end='')
            print(',%g'%(flux[i]))           
        else:
            print(',',flux[i])
          
            
def conserved_moieties(N, m_lab):
    # function to determine redundant rows
    # calculate nullspace
    NS = Matrix(N.T).nullspace()
    NS = np.array(NS).T
    # print(NS)
    rem = np.array( [] )  # container for rows that can be eliminated

    if NS.shape[0] != 0: # if nullspace is NOT empty, there are conserved moieties
        print('Conserved moieties present')
        # solutions are VECTORS (columns) in the nullspace
        # loop over all possible solutions (columns) and print out:

        for i in range(NS.shape[1]): # loop through null space basis vectors
            lab_m = np.where(NS[:,i] != 0) # find all non-zero entries
            lab_m = lab_m[0]

            for ii in range(len(lab_m)): # loop over all entries in the nullspace
                #print('+ (%g) %s ' % (NS[lab_m[ii], i], m_lab[lab_m[ii]]), end='')
                print(NS[lab_m[ii], i], m_lab[lab_m[ii]])#oben Error:Only size1-arrays can be converted to python scalars 
            print('= 0')
        
        # look which ones can go out:
        rN = linalg.matrix_rank(N)
        
        # find all possible redundant rows:
        red = np.where(np.absolute(NS.T) > 0)
        red = red[1]
        red = np.unique(red)
        # probe each row - take out and check rank
        for i in range(len(red)): # loop over the amount of dependent lines
            # if a dependent line is removed, the rank is preserved
            # check if taking out the first non-zero NS entry will work:
            # lines left after taking out redundant lines:
            rows_left = np.setdiff1d(np.arange(N.shape[0]), np.append(rem, red[i]))
            if linalg.matrix_rank(N[rows_left,:]) == rN: # check rank without balance number i
                # if rank is preserved --> removed a redundant line, store in rem
                rem = np.append(rem, red[i])
                print('Remove line for %s' % (m_lab[red[i]]))

    rem = rem.astype(int)
	
    return rem

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

def internal_cycles(N, v_lab):
    # Check for internal cycles, which cannot be determined from extracellular measurements 
    # Approach: Use nullspace analysis to find
    # all solution vectors, add constraints representing the in- and outputs of the network
    print('-'*12,'Internal cycles','-'*12)
    # Constrain extracellular fluxes
    # extracellular fluxes only have (one) input or output, not both at the same time
    n_ext_flx = np.where(np.sum(np.absolute(N), axis=0) < 1.1)
    n_ext_flx = n_ext_flx[0]

    # construct a eye matrix, select the extracellular lines:
    C_tmp = np.eye(N.shape[1])
    C = np.copy(C_tmp[n_ext_flx,:])

    NS = Matrix(np.concatenate((N,C), axis=0)).nullspace()
    NS = np.array(NS).T
    n_undef = []
    # are there fluxes in the nullspace?
    for i in range(NS.shape[0]):
        # flux in the nullspace?
        if np.any(np.absolute(NS[i,:]) > 0):
            print('flux %s cannot be determined' % (v_lab[i]))
            n_undef = np.append(n_undef, i)
    
    if len(n_undef) == 0:
        print('There are no internal cycles in this metabolic network')
    
    return n_undef

def lagrange_solve_w( E, M, Rm, std_dev_Rm ):

    # inputs:
    # matrix E (n x m)  - linearity constraints (m balance equations E . R = 0)
    # matrix M (k x m)  - measurement matrix for the k measurements (M . R = Rm +/-)
    # vector Rm (m)     - measurement values
    # vector std_dev_Rm - standard deviations of the measurement in Rm

    # rescale equation according to std_dev_Rm
    M_w = np.diag( std_dev_Rm ** -1 ) @ M
    L1 = np.hstack( (M_w.T @ M_w, E.T) )
    L2 = np.hstack( (E, np.zeros( (E.shape[0], E.shape[0]) ) ) )
    L = np.vstack( (L1, L2) )
    b = np.concatenate( (M_w.T @ (Rm/std_dev_Rm), np.zeros((E.shape[0]))))

    Rl = linalg.solve(L, b)
    R = Rl[0:E.shape[1]]

    # error propagation
    M = np.pad( M_w, ((0,0),(0,b.shape[0]-M_w.shape[1])), mode='constant')
    J = linalg.inv(L) @ M.T
    S_R = J @ J.T
    std_dev_R = np.diag(S_R) ** 0.5
    std_dev_R = std_dev_R[0:E.shape[1]]

    return (R, std_dev_R)
	
def lagrange_solve( E, M, Rm ):

    # inputs:
    # matrix E (n x m)  - linearity constraints (m balance equations E . R = 0)
    # matrix M (k x m)  - measurement matrix for the k measurements (M . R = Rm +/-)
    # vector Rm (m)     - measurement values

    L1 = np.hstack( (M.T @ M, E.T) )
    L2 = np.hstack( (E, np.zeros( (E.shape[0], E.shape[0]) ) ) )
    L = np.vstack( (L1, L2) )
    b = np.concatenate( (M.T @ Rm, np.zeros((E.shape[0]))))

    Rl = linalg.solve(L, b)
    R = Rl[0:E.shape[1]]

    return R