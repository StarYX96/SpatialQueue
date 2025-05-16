#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import random
import time
from scipy import sparse
from scipy.special import comb
from scipy.sparse.linalg import spsolve
# from dispatch_main import *
# from HomogeneousRate import *
import matplotlib.pyplot as plt

def Random_Fraction(J, seed=9001):
    '''
    :param J: Number of geographical atoms
    :return: fraction of region-wide demand f
    '''
    np.random.seed(seed)
    f = np.random.random(size=J)
    f = f / sum(f)
    return f

def Random_Pref(N, J, seed=9001):
    '''
    :param N: Number of units
    :return: Random preference list of atoms
    '''
    np.random.seed(seed=seed)
    demandNodes = np.random.rand(J, 2)
    np.random.seed(seed=seed+1)
    servers = np.random.rand(N, 2)
    PreList = np.zeros((J, N))
    for i in range(J):
        distance = np.abs(demandNodes[i] - servers).sum(axis=1)
        PreList[i] = np.argsort(distance)
    return PreList.astype(int)

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

def Myopic_Policy(N, J, Pre_L):
    # Obtain the policy to dispatch the cloest available unit given the time matrix
    policy = np.zeros([2 ** N, J], dtype=int)
    for s in range(2 ** N):
        for j in range(J):
            pre = Pre_L[j]
            for n in range(N):
                if not s >> pre[n] & 1:  # n th choice is free
                    policy[s, j] = pre[n]
                    break
    return policy

def SimulatorExp(N, J, Lambda, Mu, f, policy, uptill):
    # list indicating call locations
    k_list = np.array(random.choices(population=list(range(J)), weights=f, k=uptill))
    # list of busy units
    busy_list = []
    # Num lost call, and over threshold
    num_lost = 0
    # Steady State Distribtion
    steady_state = np.zeros(2 ** N)
    # Initialization
    inter_arrival = np.random.exponential(1 / Lambda, size=uptill)  # The inter arrival times for all arrivals
    arrival = np.cumsum(inter_arrival)  # The arrival times for all arrivals
    service_time = []
    # List of completion times for each call
    completion = np.array([[1e5 * uptill, -1]])
    # List of service completion order
    arrival_cursor = 0  # cursor in the arrival list
    completion_cursor = 0  # cursor in the completion list
    ################# Start simulation ################
    time = 0  # Starting at time 0
    while (completion_cursor + num_lost <= uptill - 1):
        while (arrival_cursor <= uptill - 1):
            ## There is an arrival at this time
            # first add the time to the past state
            state = sum(2 ** np.array(busy_list))  # current state number
            # Get preference list
            k = k_list[arrival_cursor]  # This shows the location of the call
            dis_unit = policy[state, k]  # the unit to be dispatched
            if arrival[arrival_cursor] < completion[:, 0].min():
                steady_state[state] += arrival[arrival_cursor] - time  # Update the time to stay at this state
                # Update the time to this arrival
                time = arrival[arrival_cursor]
                if state < 2 ** N - 1:  # if not last state
                    service_time.append(np.random.exponential(1 / Mu[dis_unit]))
                    # random.expovariate(1 / Mu[dis_unit])
                    completion = np.append(completion, [[arrival[arrival_cursor] + service_time[-1], dis_unit]], axis=0)
                    busy_list += [dis_unit]
                else:
                    num_lost += 1
                # Go to next arrival
                arrival_cursor += 1  # Go to next arrival
            else:
                break
        #### If this next arrival exceeds the previous , then we move to service completion
        ### Service completion
        state = sum(2 ** np.array(busy_list))  # current state number
        service_complete = completion[:, 0].argmin()
        steady_state[state] += completion[service_complete, 0] - time
        # steady_state = time_add(steady_state, busy_list, completion_time[completion_order[completion_cursor]]-time)
        # Update the time to this arrival
        time = completion[service_complete, 0]
        # If the previous call was served, take the unit out of the busy list, else pass
        try:
            busy_list.remove(completion[service_complete, 1])
        except:  # For those require 2 units
            pass
        # print('busy_list:', busy_list)
        completion = np.delete(completion, service_complete, axis=0)
        completion_cursor += 1
    # print(serve_list)
    # print(arrival)
    prob_dist = steady_state / sum(steady_state)

    return prob_dist

if __name__ == '__main__':
    Data = {}
    Data['J'] = 71
    J = Data['J']
    convCritBnD = 1e-3
    uptill = int(1e5)
    repNum = 100

    Table = pd.DataFrame(columns=['N', 'rho', 'dist_err', 'MRT_exp'])
    Table_dist = pd.DataFrame(columns=['N', 'rho', 'simu_dist_exp', 'BnD_dist'])

    for N in range(8, 13):
        for rho in [0.1, 0.5, 0.9]:
            err = []
            print('N:', N, 'rho', rho)
            np.random.seed(seed=4)
            Data['Mu'] = np.random.uniform(-0.2, 0.2, N) + 1.76
            Mu = Data['Mu']
            Data['Lambda'] = rho * Mu.sum()
            Lambda = Data['Lambda']
            Pre_L = Random_Pref(N, J, seed=4)
            Data['f'] = Random_Fraction(J, seed=4)
            f = Data['f']

            transup, transdown = Transition_Rate(N, J, Lambda, Mu)
            prob_dist_lp, calcuTime_lp = Hypercube_lp(N, transup, transdown)

            # # Simulation Solution
            policy = Myopic_Policy(N, J, Pre_L)
            simu_exp_dist = np.zeros((repNum, 2**N))

            # simuError = np.zeros(repNum)
            start_time = time.time()
            for i in range(repNum):
                simu_dist = SimulatorExp(N, J, Lambda, Mu, f, policy, uptill * 2 ** (N-9))
                simu_exp_dist[i] = simu_dist
                err.append(max(np.abs(simu_dist - prob_dist_lp) * 100 / prob_dist_lp))

            calcuTime_simu = time.time() - start_time
            print("------Simulation Time: %s seconds ------" % calcuTime_simu)
            simu_dist_avg = simu_exp_dist.mean(axis=0)
            err_avg = max(np.abs(simu_dist_avg - prob_dist_lp) * 100 / prob_dist_lp)
            print('------Simulation ERR AVG: %s ------' % err_avg)

            Table.loc[len(Table.index)] = [N, rho,  [np.array(err).mean(), np.array(err).std()]]
            del transup, transdown, prob_dist_BnD, simu_exp_dist

    Table.to_csv('Table-SimuCompareResults.csv', index=None)
