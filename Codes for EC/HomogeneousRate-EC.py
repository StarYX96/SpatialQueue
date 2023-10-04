# -*- coding: utf-8 -*-
'''
This Code is designed for comparing the performance between hypercube and BnD model.
Output: Table 1 in the paper
'''

import numpy as np
import pandas as pd
import random
import math
import time
from scipy.special import comb
from scipy import sparse


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
    servers = np.random.rand(N, 2)
    PreList = np.zeros((J, N))
    distanceList = np.zeros((J, N))
    for i in range(J):
        distance = np.abs(demandNodes[i] - servers).sum(axis=1)
        distanceList[i] = distance
        PreList[i] = np.argsort(distance)
    return PreList.astype(int), np.insert(distanceList, N, 10000, axis=1)

def Transition_Rate(N, J, Lambda, Mu):
    Pre_L, distance = Random_Pref(N, J, seed=4)
    f = Random_Fraction(J, seed=4)
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
        down_rate[down_startLoc: down_startLoc + (N - n)] = Mu
        down_startLoc += N - n


    down_row[down_startLoc: down_startLoc + N] = 0 | (1 << np.array(units))
    down_col[down_startLoc: down_startLoc + N] = 0
    down_rate[down_startLoc: down_startLoc + N] = Mu

    up_col[up_startLoc: up_startLoc + N] = stateNum - 1
    up_row[up_startLoc: up_startLoc + N] = (stateNum - 1) & ~(1 << np.array(units))
    up_rate[up_startLoc: up_startLoc + N] = [1] * N

    transup_rate = sparse.csc_matrix((up_rate * Lambda, (up_row, up_col)), shape=(stateNum, stateNum))
    transdown_rate = sparse.csc_matrix((down_rate, (down_row, down_col)), shape=(stateNum, stateNum))

    return transup_rate, transdown_rate

def BnD(N, Mu, Lambda, transup, transdown):
    """
    :param transup, transdown: upward and downward transition rate matrix in sparse form
    :return: probability distribution solved by BnD model with Guassian-Seidel iteration
    """

    # steady state probability of busy units
    rho = Lambda / Mu
    p_n = np.array([(rho ** j) / math.factorial(j) for j in range(N+1)])
    p_n = p_n / sum(p_n)

    # initial states probabilities
    p_n_B_new = np.zeros((2, 2 ** N))
    for i in range(2 ** N):
        w = bin(i).count('1')
        p_n_B_new[:, i] = [w, 1 / comb(N, w)]

    # start iteration
    transup = transup.tocsr()
    transdown = transdown.tocsr()
    layer = [np.where(p_n_B_new[0,:] == i)[0] for i in range(N+1)]
    errBnD_Linf = []
    ite = 0
    gap = 1
    time_list = []
    start_time = time.time()
    while gap >= 1e-3:
        p_n_B = np.copy(p_n_B_new[1,:])
        start_time_ite = time.time()
        for n in range(1, N):
            p_n_B_new[1, layer[n]] = (p_n_B_new[1, layer[n-1]] * transup[layer[n-1]][:, layer[n]] * n / rho + p_n_B_new[1,layer[n+1]] * transdown[layer[n+1]][:, layer[n]] * rho / (n + 1)) / (n * Mu + Lambda)
            time_list += [time.time() - start_time_ite]
            start_time_ite = time.time()
        ite += 1
        gap = max(np.abs(p_n_B_new[1,:] - p_n_B) / p_n_B_new[1,:])
        errBnD_Linf.append(gap)
    calcuTime = time.time() - start_time
    print("------ %s seconds, %s iterations ------" % (time.time() - start_time, ite))

    for i in range(N+1):
        p_n_B_new[1, layer[i]] =  p_n[i] * p_n_B_new[1, layer[i]]

    return p_n_B_new[1,:], ite, calcuTime

def LarsonOrigin(N, J, Pre_L, Distance, f, Lambda, Mu):
    '''
    :param Pre_L: preference list
    :param f: fraction of regional wide demand
    :return: upward and downward transition rate matrix
    '''
    stateNum = 2 ** N

    start_time = time.time()
    # Generate the tour sequence of states
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

    # Determine the weights of states before
    w = np.zeros(stateNum, dtype=int)
    for k in range(1, stateNum):
        if S[k] > S[k-1]:
            w[S[k]] = w[S[k - 1]] + 1
        else:
            w[S[k]] = w[S[k - 1]] - 1

    # storing the state transition matrix
    MAP = np.zeros(2**N + 1).astype(int)
    MAP[1] = 1
    for i in range(2, 2**N + 1):
        MAP[i] = MAP[i-1] + w[i-1] + 1
    trans_rate = np.zeros(MAP[-1])

    # Calculate the upward transition rate of each state
    for j in range(J):
        pre = Pre_L[j]
        dist = Distance[j]
        opt = pre[0]
        for i in range(0, len(S)):
            if S[i] == stateNum - 1:
                opt = N
                continue
            n0 = int(math.log(S[i] ^ S[i - 1], 2))
            if S[i] < S[i - 1] and dist[n0] < dist[opt]:
                opt = n0
            elif S[i] > S[i - 1] and n0 == opt:
                free = [j for j in range(N) if S[i] & (1 << j) == 0]
                r1 = np.where(np.isin(pre, free))[0][0]
                opt = pre[r1]
            dest = S[i] | (1 << opt)
            trans_rate[MAP[dest] + w[dest] - busy_set[dest].index(opt)] += f[j] * Lambda

    for i in range(len(MAP) - 2):
        trans_rate[MAP[i]] = Lambda + w[i] * Mu
    trans_rate[MAP[-2]] = w[-1] * Mu

    trans_calcuTime = time.time() - start_time
    print('Coefficients Rate Generation Time:', trans_calcuTime)

    # Iterate to solve the balance equation

    # steady state probability of busy units
    rho = Lambda / Mu
    p_n = np.array([(rho ** j) / math.factorial(j) for j in range(N+1)])
    p_n = p_n / sum(p_n)

    # initial states probabilities
    p_B = np.zeros(2 ** N)
    p_B[0] = p_n[0]
    for i in range(1, 2 ** N, 2):
        p_B[S[i]] = p_n[w[S[i]]] / comb(N, w[S[i]])

    # start iteration
    err_Linf = []
    ite = 0
    gap = 1
    start_time = time.time()
    while gap >= 1e-3:
        p_B_old = np.copy(p_B)
        for hp in range(2, 0, -1):
            for i in range(hp, 2**N, 2):
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
                p_B[S[i]] = relative / trans_rate[MAP[S[i]]]

        ite += 1
        gap = max(np.abs(p_B_old - p_B) / p_B)
        err_Linf.append(gap)
    prob_calcuTime = time.time() - start_time
    print("------ %s seconds, %s iterations ------" % (prob_calcuTime, ite))

    return p_B, ite, trans_calcuTime, prob_calcuTime

if __name__ == '__main__':
    seed = 4
    Data = {}
    Data['J'] = 99
    J = Data['J']
    Data['Mu'] = 1.58
    Mu = Data['Mu']
    Data['f'] = Random_Fraction(J, seed=seed)
    f = Data['f']

    Table = pd.DataFrame(
        columns=['N', 'rho', 'Tour_Time', 'OH_calcuTime', 'OH_ite', 'transMatrix_calcuTime',
                 'BnD_calcuTime', 'BnD_ite', 'BnD-HC-error'])

    for n in range(11, 21):
        for rho in np.arange(0.1, 1.0, 0.1):
            Data['rho'] = rho
            Data['N'] = n
            N = Data['N']
            Data['Lambda'] = rho * Mu * n
            Lambda = Data['Lambda']
            Data['Pre_L'], Data['Distance'] = Random_Pref(N, J, seed=seed)
            Pre_L, Distance = Data['Pre_L'], Data['Distance']
            print('Number of units: %s, utilization: %s' % (N, rho))

            print('------Original Hypercube Solution------')
            prob_dist_OH, ite_OH, TourTime, calcuTime_OH = LarsonOrigin(N, J, Pre_L, Distance, f, Lambda, Mu)

            print('------Birth and Death Solution------')
            start_time = time.time()
            transup, transdown = Transition_Rate(N, J, Lambda, Mu)
            trans_calcuTime = time.time() - start_time
            print('Coefficients Rate Generation Time:', trans_calcuTime)
            prob_dist_BnD, ite_BnD, calcuTime_BnD = BnD(N, Mu, Lambda, transup, transdown)

            BnDOHError = max(np.abs(prob_dist_BnD - prob_dist_OH) / prob_dist_OH)
            print('Error rate of BnD-HC model:', BnDOHError)

            Table.loc[len(Table.index)] = [N, rho, TourTime, calcuTime_OH, ite_OH, trans_calcuTime,
                                             calcuTime_BnD, ite_BnD, BnDOHError]

            del transup, transdown, prob_dist_OH, prob_dist_BnD

        Table.to_csv('Table-BnD-OHComparison.csv', index=None)