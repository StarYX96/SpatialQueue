# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import random
import math
import time
from scipy import sparse
from scipy.special import comb
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

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

# Function to generate a random preference list of atoms
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

    # Initialize arrays to store preference list
    PreList = np.zeros((J, N))

    # Loop over geographical atoms
    for i in range(J):
        # Calculate the Manhattan distances between the current geographical atom and all servers
        distance = np.abs(demandNodes[i] - servers).sum(axis=1)

        # Generate the preference list of the servers based on the distances from near to far
        PreList[i] = np.argsort(distance)
    
    return PreList

# Function to generate transition rates by traversing the states
def Transition_Rate(N, J, Lambda, Mu):
    '''
    :param N: Number of units
    :param J: Number of geographical atoms
    :param Lambda: overall arrival rate within the entire region
    :param Mu: service rate of each response unit 
    :return: Random preference list of atoms
    '''
    seed = 4
    # Generate a random preference list
    Pre_L = Random_Pref(N, J, seed=seed)

    # Generate a random fraction of region-wide demand
    f = Random_Fraction(J, seed=seed)
    
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

        # Populate column i of the upward transition matrix
        up_col[up_startLoc: up_startLoc + n] = i
        
        # Find all unit-step reachable states by changing exactly one unit from 1 to 0 of the current state
        up_row[up_startLoc: up_startLoc + n] = i & ~(1 << np.array(busy))

        # Loop over geographical atoms
        for j in range(J):

            # For each atom, find the optimal dispatch unit and add the fraction of region-wide demand to the corresponding location in upward transition matrix
            for m in range(n):
                if Pre_L[j][m] in free:
                    break
                else:
                    up_rate[up_startLoc + busy.index(Pre_L[j][m])] += f[j]
        up_startLoc += n

        # Populate column i of the downward transition matrix
        # Find all unit-step reachable states by changing exactly one unit from 0 to 1 of the current state
        down_row[down_startLoc: down_startLoc + (N - n)] = i | (1 << np.array(list(free)))
        down_col[down_startLoc: down_startLoc + (N - n)] = i

        # Fill corresponding Mu[i] in the downward transition matrix
        for m in free:
            down_rate[down_startLoc] = Mu[m]
            down_startLoc += 1

    # Populate the first column of downward transitions 
    down_row[down_startLoc: down_startLoc + N] = 0 | (1 << np.array(units))
    down_col[down_startLoc: down_startLoc + N] = 0
    for m in units:
        down_rate[down_startLoc] = Mu[m]
        down_startLoc += 1

    # Populate the last column of upward transition matrices
    up_col[up_startLoc: up_startLoc + N] = stateNum - 1
    up_row[up_startLoc: up_startLoc + N] = (stateNum - 1) & ~(1 << np.array(units))
    up_rate[up_startLoc: up_startLoc + N] = [1] * N

    # Create sparse matrices for upward and downward transition rates
    transup_rate = sparse.coo_matrix((up_rate * Lambda, (up_row, up_col)), shape=(stateNum, stateNum))
    transdown_rate = sparse.coo_matrix((down_rate, (down_row, down_col)), shape=(stateNum, stateNum))

    return transup_rate, transdown_rate

# Function to solve linear programming problem for the hypercube
def Hypercube_lp(N, transup, transdown):
    '''
    :param N: number of units
    :param transup, transdown: upward and downward transition rate matrix in sparse form
    :return: probability distribution solved by the sparse solver

    Note: This function can only solve problems with N less than 16
    '''

    # Calculate the total inflow for each state and the outflow for each state
    inflow = (transup + transdown).T
    outflow = np.array(inflow.sum(axis=0))[0]

    # Create a diagonal matrix with outflow values
    outflow_diag = sparse.diags(outflow, 0)

    # Calculate the transition matrix by subtracting outflow from inflow
    transition = inflow - outflow_diag
    transition = transition.tolil()

    # Set the last row of the transition matrix to represent the constraint that conditional probabilities sum to 1
    transition[-1] = np.ones(2 ** N)

    # Set the last element of the right-hand side to 1 to represent the sum-to-1 constraint
    b = np.zeros(2 ** N)
    b[-1] = 1

    # Solve the linear programming problem using python sparse solver
    start_time = time.time()

    # Convert to compressed sparse column format for efficient solving
    transition = transition.tocsc()
    prob_dist = spsolve(transition, b)
    calcuTime = time.time() - start_time
    print("------Solver Solution Time: %s seconds ------" % calcuTime)

    return prob_dist, calcuTime

# Function to solve the BnD model
def BnD(N, Mu, Lambda, transup, transdown, convCritBnD):
    """
    :param N: Number of units
    :param Lambda: overall arrival rate within the entire region
    :param Mu: service rate of each response unit
    :param transup, transdown: upward and downward transition rate matrix in sparse form
    :param convCritBnD: stopping criterion 
    :return: probability distribution solved by BnD model
    """
    
    # Convert transition rate matrices to CSR format for efficient computations
    transup = transup.tocsr()
    transdown = transdown.tocsr()

    # Initialize states probabilities by assuming all states within a layer is equally distributed
    p_n_B_new = np.zeros((2, 2 ** N))
    for i in range(2 ** N):
        w = bin(i).count('1')
        p_n_B_new[:, i] = [w, 1 / comb(N, w)]
        # p_B_new[:, i] = [w, 1 / (2 ** N)]

    # Calculate the downward transition rate of each state
    downRate = transdown.sum(axis=1).A.reshape(2 ** N, )

    # Define layers based on the number of busy units
    layer = [np.where(p_n_B_new[0, :] == i)[0] for i in range(N + 1)]
    
    # Initialize the downward transition rate of each layer
    MuN = np.zeros(len(Mu) + 1)
    for n in range(1, N + 1):
        MuN[n] = downRate[layer[n]].sum() / comb(N, n)

    # Initialize variables for iteration and convergence check
    errBnD_Linf = []
    ite = 0
    gap = 1
    time_list = []
    start_time = time.time()

    # Perform iterations until convergence
    while gap >= convCritBnD:
        p_n_B = np.copy(p_n_B_new[1, :])
        start_time_ite = time.time()

        # Iterate over layers using the updated function
        for n in range(1, N):
            # Repeatedly update each layer for a constant time to achieve stability
            for r in range(10):
                p_n_B_new[1, layer[n]] = (p_n_B_new[1, layer[n - 1]] * transup[layer[n - 1]][:, layer[n]] * MuN[
                    n] / Lambda + p_n_B_new[1, layer[n + 1]] * transdown[layer[n + 1]][:, layer[n]] * Lambda / MuN[
                                              n + 1]) / (downRate[layer[n]] + Lambda)
                MuN[n] = np.dot(p_n_B_new[1, layer[n]], downRate[layer[n]])

            time_list += [time.time() - start_time_ite]
            start_time_ite = time.time()

        # Update iteration count and calculate the maximum relative change
        ite += 1
        gap = max(np.abs(p_n_B_new[1, :] - p_n_B) / p_n_B_new[1, :])
        errBnD_Linf.append(gap)

    # Record the total calculation time
    calcuTime = time.time() - start_time
    print("------Birth and Death Solution Time: %s seconds, %s iterations ------" % (time.time() - start_time, ite))

    # Calculate the steady state probabilities of different layers
    p_n = np.ones(N + 1)
    for i in range(1, N + 1):
        p_n[i] = (Lambda ** i) / np.prod(MuN[1:i + 1])
    p_n = p_n / sum(p_n)

    # Recover the probability distribution
    for i in range(N + 1):
        p_n_B_new[1, layer[i]] = p_n[i] * p_n_B_new[1, layer[i]]

    return p_n_B_new[1, :], ite, calcuTime


if __name__ == '__main__':
    # Set a random seed for reproducibility
    seed = 4

    # Initialize a dictionary to store critical data
    Data = {}

    # Set the number of geographical atoms
    Data['J'] = 71
    J = Data['J']

    # Set the convergence criterion
    convCritBnD = 1e-3
    
    # Initialize a DataFrame to store results
    Table = pd.DataFrame(columns=['N', 'rho', 'LP_calcuTime', 'BnD_calcuTime', 'BnD_ite', 'BnD-error'])

    # Loop over different numbers of units and utilization rates
    for n in range(11, 16):
        for rho in np.arange(0.2, 1.0, 0.1):
            # Set data for the current iteration
            Data['rho'] = rho
            Data['N'] = n
            N = Data['N']

            # Set the heterogeneous service rates
            np.random.seed(seed=seed)
            Data['Mu'] = np.random.uniform(-0.2, 0.2, n) + 1.76
            Mu = Data['Mu']

            # Set overall arrival rate within entire region
            Data['Lambda'] = rho * Mu.sum()
            Lambda = Data['Lambda']
            
            print('Number of units: %s, utilization: %s' % (N, rho))
            
            # Generate transition rates
            transup, transdown = Transition_Rate(N, J, Lambda, Mu)
            
            # Solve the Linear Programming Model 
            prob_dist_lp, calcuTime_lp = Hypercube_lp(N, transup, transdown)

            # Solve the Birth and Death Model
            prob_dist_BnD, ite_BnD, calcuTime_BnD = BnD(N, Mu, Lambda, transup, transdown, convCritBnD)

            # Calculate the error rate between two models
            BnDError = max(np.abs(prob_dist_BnD - prob_dist_lp) / prob_dist_lp)
            print('Relative Error:', BnDError)

            # Record results in the DataFrame
            Table.loc[len(Table.index)] = [N, rho, calcuTime_lp, calcuTime_BnD, ite_BnD, BnDError]

            # Clean up variables and release the memory for next iteration
            del transup, transdown, prob_dist_lp, prob_dist_BnD

    # Save the results to a CSV file
    Table.to_csv('Table-BnD-SolverComparison.csv', index=None)



