# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import random
import math
import pickle
import time
import itertools
from scipy import sparse
from scipy.special import comb
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt


def Transition_Rate(N, J, Lambda, Mu, Pre_L, f):
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
    # initial states probabilities
    p_n_B_new = np.zeros((2, 2 ** N))
    for i in range(2 ** N):
        w = bin(i).count('1')
        p_n_B_new[:, i] = [w, 1 / comb(N, w)]

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
    p_n_B_new[1] = prob_dist
    calcuTime = time.time() - start_time
    print("------Solver: %s seconds ------" % calcuTime)

    return p_n_B_new, calcuTime

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
    # print("------Birth and Death Solution Time: %s seconds, %s iterations ------" % (time.time() - start_time, ite))

    p_n = np.ones(N + 1)
    for i in range(1, N + 1):
        p_n[i] = (Lambda ** i) / np.prod(MuN[1:i + 1])
    p_n = p_n / sum(p_n)
    for i in range(N + 1):
        p_n_B_new[1, layer[i]] = p_n[i] * p_n_B_new[1, layer[i]]

    return p_n_B_new, ite, calcuTime

def MRT(N, J, prob, pre, f, t, threshold):
    rho = np.zeros((J, N))
    frac = np.zeros(J + 1)
    for i in range(2 ** N - 1):
        busy = [m for m in range(N) if i & (1 << m) != 0]
        for j in range(J):
            n = int(next((m for m in pre[j] if m not in busy), None))
            rho[j, n] += prob[1, i] * f[j]
            if t[j, n] <= threshold:
                frac[j] += prob[1, i]
                frac[-1] += prob[1, i] * f[j]
    rho = rho / (1 - prob[1, -1])
    frac = frac / (1 - prob[1, -1])
    MRT = t * rho

    return MRT.sum(), frac


def parameter(File_Name='0_2', Lambda=None, Mu=None, e=None, f=None, Pre_L_1=None, Pre_L_2=None, sort=True):
    '''
		:File_Name: File that stored the parameter for the runs
		:Lambda, Mu, e, f, Pre_L_1, Pre_L_2: the date of this collection
		Output: Data
	'''
    file = open(File_Name + '.pkl', 'rb')
    Data = pickle.load(file)
    if Lambda != None:
        Data['Lambda_1'], Data['Lambda_2'] = Lambda
    if Mu != None:
        Data['Mu_1'], Data['Mu_2'] = Mu
    if e != None:
        Data['e'] = np.array(e)
    if f != None:
        Data['f'] = np.array(f)
    if np.any(Pre_L_1 != None):
        Data['Pre_L_1'] = Pre_L_1
    if np.any(Pre_L_2 != None):
        Data['Pre_L_2'] = Pre_L_2
    N, N_1, N_2 = Data['N'], Data['N_1'], Data['N_2']
    if sort:  # If want to perform argsort
        Data['Pre_L_1'] = Data['Pre_L_1'].argsort().argsort()
        Data['Pre_L_2'] = Data['Pre_L_2'].argsort().argsort()
    file.close()
    return Data

def Data_to_Param(Data):
    N_1, N_2 = Data['N_1'], Data['N_2']
    K = Data['K']
    Lambda_1, Lambda_2 = Data['Lambda_1'], Data['Lambda_2']
    Mu_1, Mu_2 = Data['Mu_1'], Data['Mu_2']
    N = Data['N']
    pre_list_1, pre_list_2 = Data['Pre_L_1'], Data['Pre_L_2']
    frac_j_1, frac_j_2 = Data['e'], Data['f']
    return N_1, N_2, K, Lambda_1, Lambda_2, Mu_1, Mu_2, N, pre_list_1, pre_list_2, frac_j_1, frac_j_2

if __name__ == '__main__':
    seed = 4
    Data = parameter(File_Name='Saint_Paul')
    N_1, N_2, K, Lambda_1, Lambda_2, Mu_1, Mu_2, N, pre_list_1, pre_list_2, frac_j_1, frac_j_2 = Data_to_Param(Data)
    N = 9
    J = K
    Lambda = Lambda_1
    np.random.seed(seed=seed)
    Mu = np.random.uniform(0.8, 1.2, N) * Mu_1
    convCritBnD = 1e-3

    distance = pd.read_csv('distance_EMS.csv')
    distance = distance.drop(distance.columns[0], axis=1)
    distance = distance / 1600
    distance[distance > 0.7712735431536174] = (distance[distance > 0.7712735431536174] * 111.51113889304331 + 86.005591195132666) / 60
    distance[distance < 0.7712735431536174] = 195.86302790816589 * np.sqrt(distance[distance < 0.7712735431536174]) / 60
    t = np.array(distance)
    Pre_L = np.zeros((J, N))
    for i in range(J):
        Pre_L[i] = np.argsort(t[i])
    Pre_L = Pre_L.astype(int)
    f = frac_j_1
    threshold = 4.5

    Table = pd.DataFrame(columns=['PreList', 'prob_BnD', 'utilization', 'MRT_6', 'frac_total_6', 'MRT_9', 'frac_total_9'])

    vary = 4

    PreList = [Pre_L]
    for j in range(J):
        pre_j = np.copy(Pre_L)
        JVary = itertools.permutations(Pre_L[j][:vary])
        for pre in JVary:
            pre_j[j][:vary] = np.array(pre)
            PreList.append(np.copy(pre_j))

    print('Total cases:', len(PreList))
    for pre in PreList:
        transup, transdown = Transition_Rate(N, J, Lambda, Mu, pre, f)

        prob_dist_BnD, ite_BnD, calcuTime_BnD = BnD(N, Mu, Lambda, transup, transdown, convCritBnD)

        utilizeRate = np.zeros(N)
        for i in range(2 ** N):
            n = '{0:09b}'.format(i)
            for j in range(N):
                if n[-j - 1] == '1':
                    utilizeRate[j] += prob_dist_BnD[1, i]

        MRT_BnD_6, thresholdFrac_6 = MRT(N, J, prob_dist_BnD, pre, f, t, threshold)
        MRT_BnD_9, thresholdFrac_9 = MRT(N, J, prob_dist_BnD, pre, f, t, threshold + 3)

        Table.loc[len(Table.index)] = [pre, prob_dist_BnD[1], utilizeRate, MRT_BnD_6, thresholdFrac_6[-1], MRT_BnD_9, thresholdFrac_9[-1]]

        del transup, transdown, prob_dist_BnD, utilizeRate

    Table.to_pickle('Table-OptDispatchPolicy.pkl')


