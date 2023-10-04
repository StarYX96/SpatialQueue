# -*- coding: utf-8 -*-
'''
This Code is designed for comparing the performance between hypercube and BnD model.
Output: Table 1 in the paper
'''

import numpy as np
import pandas as pd
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


def LarsonSparse(N, J, Pre_L, Distance, f, Lambda, Mu):
    '''
    :param Pre_L: preference list
    :param f: fraction of regional wide demand
    :return: upward and downward transition rate matrix
    '''
    stateNum = 2 ** N

    # Generate the tour sequence of states
    start_time = time.time()
    busy_set = [[0]] * (2 ** N)
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
    print('Transition Rate Calculation Time:', trans_calcuTime)

    # steady state probability of busy units
    rho = Lambda / Mu
    p_n = np.array([(rho ** j) / math.factorial(j) for j in range(N + 1)])
    p_n = p_n / sum(p_n)

    # initial states probabilities
    p_B = np.zeros((2, 2 ** N))
    for i in range(2 ** N):
        w = bin(i).count('1')
        p_B[:, i] = [w, p_n[w] / comb(N, w)]

    # recover the MAP to coo sparse matrix
    start_time = time.time()
    up_row = np.zeros(len(trans_rate), dtype=int)
    up_col = np.zeros(len(trans_rate), dtype=int)
    startLoc = 0
    for i in range(len(MAP) - 1):
        trans_rate[MAP[i]] = 0
        up_row[startLoc: MAP[i+1] - 1] = np.array([i & ~(1 << j) for j in busy_set[i][::-1]])
        up_col[startLoc: MAP[i+1] - 1] = i
        startLoc = MAP[i+1]

    transup = sparse.csc_matrix((trans_rate[1:], (up_row[:-1], up_col[:-1])), shape=(stateNum, stateNum))
    noise = np.where(up_col == 0)[0]
    down_row = np.delete(up_col, noise)
    down_col = np.delete(up_row, noise)
    down_rate = np.array([Mu] * len(down_row))
    transdown = sparse.csc_matrix((down_rate, (down_row, down_col)), shape=(stateNum, stateNum))

    # start iteration
    transup = transup.tocsr()
    transdown = transdown.tocsr()
    transition = transup + transdown
    odd = np.where(p_B[0,] % 2 == 1)[0]
    even = np.where(p_B[0,] % 2 == 0)[0]
    ite = 0
    gap = 1
    while gap >= 1e-3:
        p_B_old = np.copy(p_B[1, :])
        p_B[1, odd] = ((p_B[1, even] * transition[even][:, odd]) / (Lambda + p_B[0, odd] * Mu))
        p_B[1, -1] = p_B_old[-1]
        p_B[1, even[1:]] = ((p_B[1, odd] * transition[odd][:, even]) / (Lambda + p_B[0, even] * Mu))[1:]
        p_B[1, -1] = p_B_old[-1]
        ite += 1
        gap = max(np.abs(p_B_old - p_B[1, :]) / p_B[1, :])
    prob_calcuTime = time.time() - start_time
    print("------Modified Hypercube Solution: %s seconds, %s iterations ------" % (prob_calcuTime, ite))
    return p_B[1, :], trans_calcuTime, prob_calcuTime

if __name__ == '__main__':
    seed = 4
    Data = {}
    Data['J'] = 71
    J = Data['J']
    Data['Mu'] = 1.76
    Mu = Data['Mu']
    Data['f'] = Random_Fraction(J, seed=seed)
    f = Data['f']

    Table = pd.DataFrame(columns=['N', 'rho', 'transMatrix_calcuTime', 'Prob_calcuTime'])

    for n in range(11, 21):
        for rho in np.arange(0.1, 1.0, 0.1):
            Data['rho'] = rho
            Data['N'] = n
            N = Data['N']
            Data['Lambda'] = rho * Mu * n
            Lambda = Data['Lambda']
            Data['Pre_L'], Data['Distance'] = Random_Pref(N, J, seed=seed)
            Pre_L, Distance = Data['Pre_L'], Data['Distance']

            print('Number of units %s, utilization: %s' %(N, rho))
            prob_dist, trans_calcuTime, prob_calcuTime = LarsonSparse(N, J, Pre_L, Distance, f, Lambda, Mu)
            Table.loc[len(Table.index)] = [N, rho, trans_calcuTime, prob_calcuTime]
            del prob_dist

    Table.to_csv('Table-HomogeneousRate-ModifiedHypercube.csv', index=None)