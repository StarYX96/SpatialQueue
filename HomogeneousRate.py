# -*- coding: utf-8 -*-

'''
This code is designed for evaluating the performance of oringal hypercube model and birth anf death model.
'''

import numpy as np
import pandas as pd
import random
import math
import time
from scipy.special import comb
from scipy import sparse

# Function to generate a random fraction of region-wide demand
def Random_Fraction(J, seed=9001):
    '''
    :param J: Number of geographical atoms
    :return: fraction of region-wide demand f
    '''
    # Set the random seed 
    np.random.seed(seed)
    
    # Generate an array of random numbers between 0 and 1
    f = np.random.random(size=J)

    # Normalize the array to represent a fraction (summing to 1)
    f = f / sum(f)
    return f

# Function to generate a random preference list of atoms and their distances
def Random_Pref(N, J, seed=9001):
    '''
    :param N: Number of units
    :param J: Number of geographical atoms
    :return: Random preference list of atoms
    '''
    np.random.seed(seed=seed)

    # Generate random coordinates for geographical atoms
    demandNodes = np.random.rand(J, 2)

    # Generate random coordinates for servers
    servers = np.random.rand(N, 2)

    # Initialize arrays to store preference list and distances of atom-server pairs
    PreList = np.zeros((J, N))
    distanceList = np.zeros((J, N))

    # Loop over geographical atoms
    for i in range(J):
        # Calculate and store the Manhattan distances between the current geographical atom and all servers
        distance = np.abs(demandNodes[i] - servers).sum(axis=1)
        distanceList[i] = distance
        
        # Generate the preference list of the servers based on the distances from near to far
        PreList[i] = np.argsort(distance)

    # Convert the preference list to integer type and construct an artificial server with a large distance (10000) for computational convenience
    return PreList.astype(int), np.insert(distanceList, N, 10000, axis=1)

# Function to generate transition rates by traversing the states  
def Transition_Rate(N, J, Lambda, Mu):
    '''
    :param N: Number of units
    :param J: Number of geographical atoms
    :param Lambda: overall arrival rate within entire region
    :param Mu: service rate of each response unit 
    :return: Random preference list of atoms
    '''
    
    # Generate a random preference list and distances
    Pre_L, distance = Random_Pref(N, J, seed=4)
    # Generate a random fraction of region-wide demand
    f = Random_Fraction(J, seed=4)
    
    stateNum = 2 ** N
    units = list(range(N))

    # Calculate the density of the transition matrix
    density = N * (2 ** (N - 1))

    # Initialize arrays to store the upward and downward transition rates
    up_row = np.zeros(density)
    up_col = np.zeros(density)
    up_rate = np.zeros(density)
    down_row = np.zeros(density)
    down_col = np.zeros(density)
    down_rate = np.zeros(density)

    # Initialize counters for the starting location of upward and downward transitions in each loop
    up_startLoc = 0
    down_startLoc = 0

    # Loop over all states (except the first and last)
    for i in range(1, stateNum - 1):
        # Identify busy units in the current state
        busy = [m for m in range(N) if i & (1 << m) != 0]
        n = len(busy)
        # Identify free units in the current state
        free = set(units) - set(busy)

        # Populate the column i of upward transition matrix
        up_col[up_startLoc: up_startLoc + n] = i

        # Find all unit-step reachable states by changing exactly one units from 1 to 0 of current state
        up_row[up_startLoc: up_startLoc + n] = i & ~(1 << np.array(busy))

        # Loop over geographical atoms
        for j in range(J):
            
            # For each atom, find the optimal dispatch unit and add the fraction of region-wide demand to corresponding location in upward transition matrix
            for m in range(n):
                if Pre_L[j][m] in free:
                    break
                else:
                    up_rate[up_startLoc + busy.index(Pre_L[j][m])] += f[j]
        up_startLoc += n

        # Populate the column i of downward transition matrix
        # Find all unit-step reachable states by changing exactly one units from 0 to 1 of current state
        down_row[down_startLoc: down_startLoc + (N - n)] = i | (1 << np.array(list(free)))
        down_col[down_startLoc: down_startLoc + (N - n)] = i

        # Fill ``Mu''  to corresponding locations in the downward transition matrix
        down_rate[down_startLoc: down_startLoc + (N - n)] = Mu
        down_startLoc += N - n


    # Populate the first column of downward transitions 
    down_row[down_startLoc: down_startLoc + N] = 0 | (1 << np.array(units))
    down_col[down_startLoc: down_startLoc + N] = 0
    down_rate[down_startLoc: down_startLoc + N] = Mu

    # Populate the last column of upward transition matrices
    up_col[up_startLoc: up_startLoc + N] = stateNum - 1
    up_row[up_startLoc: up_startLoc + N] = (stateNum - 1) & ~(1 << np.array(units))
    up_rate[up_startLoc: up_startLoc + N] = [1] * N

    # Create sparse matrices for upward and downward transition rates
    transup_rate = sparse.coo_matrix((up_rate * Lambda, (up_row, up_col)), shape=(stateNum, stateNum))
    transdown_rate = sparse.coo_matrix((down_rate, (down_row, down_col)), shape=(stateNum, stateNum))

    return transup_rate, transdown_rate

# Function to solve the birth and death model 
def BnD(N, Mu, Lambda, transup, transdown):
    """
    :param N: Number of units
    :param Lambda: overall arrival rate within entire region
    :param Mu: service rate of each response unit 
    :param transup, transdown: upward and downward transition rate matrix in sparse form
    :return: probability distribution solved by BnD model with Guassian-Seidel iteration
    """

    # Calculate the steady state probability of busy units
    rho = Lambda / Mu
    p_n = np.array([(rho ** j) / math.factorial(j) for j in range(N+1)])
    p_n = p_n / sum(p_n)

    # Initialize states probabilities by assuming all states within a layer is equally distributed
    p_n_B_new = np.zeros((2, 2 ** N))
    for i in range(2 ** N):
        w = bin(i).count('1')
        p_n_B_new[:, i] = [w, 1 / comb(N, w)]

    # Convert transition rate matrices to CSR format for efficient computations
    transup = transup.tocsr()
    transdown = transdown.tocsr()

    # Define layers based on the number of busy units
    layer = [np.where(p_n_B_new[0,:] == i)[0] for i in range(N+1)]

    # Initialize variables for iteration and convergence check
    errBnD_Linf = []
    ite = 0
    gap = 1
    time_list = []
    start_time = time.time()

    # Perform iterations until convergence
    while gap >= 1e-3:
        p_n_B = np.copy(p_n_B_new[1,:])
        start_time_ite = time.time()

        # Iterate over layers using updated function
        for n in range(1, N):
            p_n_B_new[1, layer[n]] = (p_n_B_new[1, layer[n-1]] * transup[layer[n-1]][:, layer[n]] * n / rho + p_n_B_new[1,layer[n+1]] * transdown[layer[n+1]][:, layer[n]] * rho / (n + 1)) / (n * Mu + Lambda)
            time_list += [time.time() - start_time_ite]
            start_time_ite = time.time()

        # Update iteration count and calculate the maximum relative change
        ite += 1
        gap = max(np.abs(p_n_B_new[1,:] - p_n_B) / p_n_B_new[1,:])
        errBnD_Linf.append(gap)

    # Record the total calculation time
    calcuTime = time.time() - start_time
    print("------ %s seconds, %s iterations ------" % (time.time() - start_time, ite))

    # Recover the probability distribution
    for i in range(N+1):
        p_n_B_new[1, layer[i]] =  p_n[i] * p_n_B_new[1, layer[i]]

    return p_n_B_new[1,:], ite, calcuTime

# Function to solve Larson's original hypercube model
def LarsonOrigin(N, J, Pre_L, Distance, f, Lambda, Mu):
    '''
    :param N: Number of units
    :param J: Number of geographical atoms
    :param Pre_L: preference list
    :param Distance: distance of each (atom, server) pair
    :param f: fraction of regional wide demand
    :param Lambda: overall arrival rate within entire region
    :param Mu: service rate of each response unit 
    :return: upward and downward transition rate matrix
    '''
    stateNum = 2 ** N
    start_time = time.time()
    
    # Generate the tour sequence of states and the set of busy units for each state
    busy_set = [0] * (2 ** N)
    busy_set[1] = [0]
    S = np.zeros(stateNum, dtype=int)
    S[1] = 1
    m2 = 2
    for n in range(1, N):
        m1 = m2
        m2 = 2 * m1
        for i in range(m1, m2):
            S[i] = m1 + S[m2 - i - 1]
            busy_set[S[i]] = sorted([j for j in range(N) if S[i] & (1 << j) != 0])

    # Determine the weights of states
    w = np.zeros(stateNum, dtype=int)
    for k in range(1, stateNum):
        if S[k] > S[k-1]:
            w[S[k]] = w[S[k - 1]] + 1
        else:
            w[S[k]] = w[S[k - 1]] - 1

    # Determine the address of the on-diagonal term in transition rate array
    MAP = np.zeros(2**N + 1).astype(int)
    MAP[1] = 1
    for i in range(2, 2**N + 1):
        MAP[i] = MAP[i-1] + w[i-1] + 1

    # Creat the transition rate array
    trans_rate = np.zeros(MAP[-1])

    # Calculate the upward transition rate of each state
    # Iterate over geographical atoms
    for j in range(J):
        pre = Pre_L[j]
        dist = Distance[j]
        opt = pre[0]

        # Iterate over states according to the tour squence
        for i in range(0, len(S)):
            if S[i] == stateNum - 1:
                opt = N
                continue
            
            # Determine the optimal unit to dispatch
            n0 = int(math.log(S[i] ^ S[i - 1], 2))
            if S[i] < S[i - 1] and dist[n0] < dist[opt]:
                opt = n0
            elif S[i] > S[i - 1] and n0 == opt:
                free = [j for j in range(N) if S[i] & (1 << j) == 0]
                r1 = np.where(np.isin(pre, free))[0][0]
                opt = pre[r1]

            # Find the storage location in the transition rate array 
            dest = S[i] | (1 << opt)
            trans_rate[MAP[dest] + w[dest] - busy_set[dest].index(opt)] += f[j] * Lambda

    # Calculate the transition rates for states without external transitions
    for i in range(len(MAP) - 2):
        trans_rate[MAP[i]] = Lambda + w[i] * Mu
    trans_rate[MAP[-2]] = w[-1] * Mu

    trans_calcuTime = time.time() - start_time
    print('Coefficients Rate Generation Time:', trans_calcuTime)

    # Iterate to solve the balance equation
    # Calculate the steady state probability of busy units
    rho = Lambda / Mu
    p_n = np.array([(rho ** j) / math.factorial(j) for j in range(N+1)])
    p_n = p_n / sum(p_n)

    # Initialize states probabilities by assuming all states within a layer is equally distributed
    p_B = np.zeros(2 ** N)
    p_B[0] = p_n[0]
    for i in range(1, 2 ** N, 2):
        p_B[S[i]] = p_n[w[S[i]]] / comb(N, w[S[i]])

    # Initialize variables for iteration and convergence check
    err_Linf = []
    ite = 0
    gap = 1
    start_time = time.time()

    # Start iteration
    while gap >= 1e-3:
        p_B_old = np.copy(p_B)
        
        #  Alternates between updating probabilities for even- and odd-numbered layers
        for hp in range(2, 0, -1):
            
            # Iterate over states within even- or odd-numbered layers
            for i in range(hp, 2**N, 2):
                # Recover the transition rates in by MAP and applied in the the iteration equations 
                address = MAP[S[i]] + w[S[i]] + 1
                relative = 0
                for k in range(N):
                    x = S[i] ^ (1 << k)
                    cbit = (S[i] >> k) & 1
                    if cbit == 0:
                        relative += p_B[x] * Mu
                    else:
                        address = address - 1
                        relative += p_B[x] * trans_rate[address]

                # Update the probability of the current state
                p_B[S[i]] = relative / trans_rate[MAP[S[i]]]

        # Update iteration count and calculate the maximum relative change
        ite += 1
        gap = max(np.abs(p_B_old - p_B) / p_B)
        err_Linf.append(gap)

    # Record the time taken for probability calculation
    prob_calcuTime = time.time() - start_time
    print("------ %s seconds, %s iterations ------" % (prob_calcuTime, ite))

    return p_B, ite, trans_calcuTime, prob_calcuTime

if __name__ == '__main__':
    # Set a random seed for reproducibility
    seed = 4

    # Initialize a dictionary to store critical data
    Data = {}

    # Set the number of geographical atoms
    Data['J'] = 71
    J = Data['J']

    # Set the service rate
    Data['Mu'] = 1.76
    Mu = Data['Mu']

    # Generate a random fraction of region-wide demand
    Data['f'] = Random_Fraction(J, seed=seed)
    f = Data['f']

    # Initialize a DataFrame to store results
    Table = pd.DataFrame(
        columns=['N', 'rho', 'Tour_Time', 'OH_calcuTime', 'OH_ite', 'transMatrix_calcuTime',
                 'BnD_calcuTime', 'BnD_ite', 'BnD-HC-error'])

    # Loop over different numbers of units and utilization rates
    for n in range(11, 21):
        for rho in np.arange(0.1, 1.0, 0.1):
            # Set data for the current iteration
            Data['rho'] = rho
            Data['N'] = n
            N = Data['N']

            # Set overall arrival rate within entire region
            Data['Lambda'] = rho * Mu * n
            Lambda = Data['Lambda']
            
            # Generate random preference list and distance list
            Data['Pre_L'], Data['Distance'] = Random_Pref(N, J, seed=seed)
            Pre_L, Distance = Data['Pre_L'], Data['Distance']

            print('Number of units: %s, utilization: %s' % (N, rho))

            # Solve the Original Hypercube Model
            print('------Original Hypercube Solution------')
            prob_dist_OH, ite_OH, TourTime, calcuTime_OH = LarsonOrigin(N, J, Pre_L, Distance, f, Lambda, Mu)

            # Solve the Birth and Death Model
            print('------Birth and Death Solution------')
            start_time = time.time()
            transup, transdown = Transition_Rate(N, J, Lambda, Mu)
            trans_calcuTime = time.time() - start_time
            print('Coefficients Rate Generation Time:', trans_calcuTime)
            prob_dist_BnD, ite_BnD, calcuTime_BnD = BnD(N, Mu, Lambda, transup, transdown)

            # Calculate the error rate between two models
            BnDOHError = max(np.abs(prob_dist_BnD - prob_dist_OH) / prob_dist_OH)
            print('Error rate of BnD-OH model:', BnDOHError)

            # Record results in the DataFrame
            Table.loc[len(Table.index)] = [N, rho, TourTime, calcuTime_OH, ite_OH, trans_calcuTime,
                                             calcuTime_BnD, ite_BnD, BnDOHError]

            # Clean up variables and release the memory for next iteration
            del transup, transdown, prob_dist_OH, prob_dist_BnD

        # Save the results to a CSV file
        Table.to_csv('Table-BnD-OHComparison.csv', index=None)
