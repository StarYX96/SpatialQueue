# -*- coding: utf-8 -*-
'''
This Code is designed for calculate the BnD model of heterogeneous service rate when rho = 0.1
r is set as follows:
11 - 13 units: r = 30;
14 - 15 units: r = 40
16 units: r = 50
17 - 19 units: r = 60
20 units: r = 70
'''

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
    np.random.seed(seed)
    f = np.random.random(size=J)
    f = f / sum(f)
    return f

# Function to generate a random preference list of atoms
def Random_Pref(N, J, seed=9001):
    '''
    :param N: Number of units
    :return: Random preference list of atoms
    '''
    np.random.seed(seed=seed)
    demandNodes = np.random.rand(J, 2)
    servers = np.random.rand(N, 2)
    PreList = np.zeros((J, N))
    for i in range(J):
        distance = np.abs(demandNodes[i] - servers).sum(axis=1)
        PreList[i] = np.argsort(distance)
    return PreList

# Function to generate transition rates
def Transition_Rate(N, J, Lambda, Mu):
    seed = 4
    Pre_L = Random_Pref(N, J, seed=seed)
    f = Random_Fraction(J, seed=seed)
    stateNum = 2 ** N
    units = list(range(N))
    density = N * (2 ** (N - 1))

    up_row = np.zeros(density)
    up_col = np.zeros(density)
    up_rate = np.zeros(density)
    down_row = np.zeros(density)
    down_col = np.zeros(density)
    down_rate = np.zeros(density)
    up_startLoc = 0
    down_startLoc = 0

    for i in range(1, stateNum - 1):
        busy = [m for m in range(N) if i & (1 << m) != 0]
        n = len(busy)
        free = set(units) - set(busy)
        up_col[up_startLoc: up_startLoc + n] = i
        up_row[up_startLoc: up_startLoc + n] = i & ~(1 << np.array(busy))
        for j in range(J):
            for m in range(n):
                if Pre_L[j][m] in free:
                    break
                else:
                    up_rate[up_startLoc + busy.index(Pre_L[j][m])] += f[j]
        up_startLoc += n

        down_row[down_startLoc: down_startLoc + (N - n)] = i | (1 << np.array(list(free)))
        down_col[down_startLoc: down_startLoc + (N - n)] = i
        for m in free:
            down_rate[down_startLoc] = Mu[m]
            down_startLoc += 1


    down_row[down_startLoc: down_startLoc + N] = 0 | (1 << np.array(units))
    down_col[down_startLoc: down_startLoc + N] = 0
    for m in units:
        down_rate[down_startLoc] = Mu[m]
        down_startLoc += 1

    up_col[up_startLoc: up_startLoc + N] = stateNum - 1
    up_row[up_startLoc: up_startLoc + N] = (stateNum - 1) & ~(1 << np.array(units))
    up_rate[up_startLoc: up_startLoc + N] = [1] * N

    transup_rate = sparse.csc_matrix((up_rate * Lambda, (up_row, up_col)), shape=(stateNum, stateNum))
    transdown_rate = sparse.csc_matrix((down_rate, (down_row, down_col)), shape=(stateNum, stateNum))

    return transup_rate, transdown_rate

# Function to solve linear programming problem for the hypercube
def Hypercube_lp(N, transup, transdown):
    '''
    :param N: number of units
    :param transup, transdown: upward and downward transition rate matrix in sparse form
    :return: probability distribution solved by linear programming

    Note: This function can only solved the problems with N less than 15
    '''
    inflow = (transup + transdown).T
    outflow = np.array(inflow.sum(axis=0))[0]
    # outflow_diag = np.diag(outflow)
    # transition = inflow.toarray() - outflow_diag
    outflow_diag = sparse.diags(outflow, 0)
    transition = inflow - outflow_diag
    transition = transition.tolil()
    transition[-1] = np.ones(2 ** N)
    b = np.zeros(2 ** N)
    b[-1] = 1

    # Solve the problem using python sparse solver
    start_time = time.time()
    # transition_sparse = sparse.csc_matrix(transition)
    transition = transition.tocsc()
    prob_dist = spsolve(transition, b)
    calcuTime = time.time() - start_time
    print("------Solver Solution Time: %s seconds ------" % calcuTime)

    return prob_dist, calcuTime

# Function to solve the BnD model
def BnD(N, Mu, Lambda, transup, transdown, convCritBnD):
    """
    :param transup, transdown: upward and downward transition rate matrix in sparse form
    :return: probability distribution solved by BnD model with Guassian-Seidel iteration
    """

    transup = transup.tocsr()
    transdown = transdown.tocsr()

    # initial states probabilities
    p_n_B_new = np.zeros((2, 2 ** N))
    for i in range(2 ** N):
        w = bin(i).count('1')
        p_n_B_new[:, i] = [w, 1 / comb(N, w)]
        # p_B_new[:, i] = [w, 1 / (2 ** N)]
    downRate = transdown.sum(axis=1).A.reshape(2 ** N, )
    MuN = np.zeros(len(Mu) + 1)

    # start iteration
    layer = [np.where(p_n_B_new[0, :] == i)[0] for i in range(N + 1)]
    for n in range(1, N + 1):
        MuN[n] = downRate[layer[n]].sum() / comb(N, n)
    errBnD_Linf = []
    ite = 0
    gap = 1
    time_list = []
    start_time = time.time()
    while gap >= convCritBnD:
        p_n_B = np.copy(p_n_B_new[1, :])
        start_time_ite = time.time()
        for n in range(1, N):
            for r in range(10):
                p_n_B_new[1, layer[n]] = (p_n_B_new[1, layer[n - 1]] * transup[layer[n - 1]][:, layer[n]] * MuN[
                    n] / Lambda + p_n_B_new[1, layer[n + 1]] * transdown[layer[n + 1]][:, layer[n]] * Lambda / MuN[
                                              n + 1]) / (downRate[layer[n]] + Lambda)
                MuN[n] = np.dot(p_n_B_new[1, layer[n]], downRate[layer[n]])

            time_list += [time.time() - start_time_ite]
            start_time_ite = time.time()

        ite += 1
        gap = max(np.abs(p_n_B_new[1, :] - p_n_B) / p_n_B_new[1, :])
        errBnD_Linf.append(gap)
    calcuTime = time.time() - start_time
    print("------Birth and Death Solution Time: %s seconds, %s iterations ------" % (time.time() - start_time, ite))

    p_n = np.ones(N + 1)
    for i in range(1, N + 1):
        p_n[i] = (Lambda ** i) / np.prod(MuN[1:i + 1])
    p_n = p_n / sum(p_n)
    for i in range(N + 1):
        p_n_B_new[1, layer[i]] = p_n[i] * p_n_B_new[1, layer[i]]

    return p_n_B_new[1, :], ite, calcuTime


if __name__ == '__main__':
    Data = {}
    Data['J'] = 71
    J = Data['J']
    convCritBnD = 1e-3
    seed = 4

    Table = pd.DataFrame(columns=['N', 'rho', 'LP_calcuTime', 'BnD_calcuTime', 'BnD_ite', 'BnD-error'])

    # Loop over different numbers of units and utilization rates
    for n in range(11, 16):
        for rho in np.arange(0.2, 1.0, 0.1):
            Data['rho'] = rho
            Data['N'] = n
            N = Data['N']
            np.random.seed(seed=seed)
            Data['Mu'] = np.random.uniform(-0.2, 0.2, n) + 1.76
            Mu = Data['Mu']
            Data['Lambda'] = rho * Mu.sum()
            Lambda = Data['Lambda']
            print('Number of units: %s, utilization: %s' % (N, rho))
            
            # Generate transition rates
            transup, transdown = Transition_Rate(N, J, Lambda, Mu)
            
            # Solve linear programming problem for the hypercube
            prob_dist_lp, calcuTime_lp = Hypercube_lp(N, transup, transdown)

            # Solve the BnD model
            prob_dist_BnD, ite_BnD, calcuTime_BnD = BnD(N, Mu, Lambda, transup, transdown, convCritBnD)

            # Calculate relative error
            BnDError = max(np.abs(prob_dist_BnD - prob_dist_lp) / prob_dist_lp)
            print('Relative Error:', BnDError)

            # Store results in a table
            Table.loc[len(Table.index)] = [N, rho, calcuTime_lp, calcuTime_BnD, ite_BnD, BnDError]

            # Clean up
            del transup, transdown, prob_dist_lp, prob_dist_BnD

    # Save the results to a CSV file
    Table.to_csv('Table-BnD-SolverComparison.csv', index=None)



