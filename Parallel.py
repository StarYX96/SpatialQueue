# -*- coding: utf-8 -*-
'''
test for GPU parallel: transition matrix
'''
import numpy as np
import pandas as pd
import random
import math
import time
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.special import comb
from scipy.sparse.linalg import spsolve
from functools import partial
import multiprocessing as mp

from multiprocessing.reduction import ForkingPickler, AbstractReducer


class ForkingPickler4(ForkingPickler):
    def __init__(self, *args):
        if len(args) > 1:
            args[1] = 2
        else:
            args.append(2)
        super().__init__(*args)

    @classmethod
    def dumps(cls, obj, protocol=4):
        return ForkingPickler.dumps(obj, protocol)


def dump(obj, file, protocol=4):
    ForkingPickler4(file, protocol).dump(obj)


class Pickle4Reducer(AbstractReducer):
    ForkingPickler = ForkingPickler4
    register = ForkingPickler4.register
    dump = dump


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
    for i in range(J):
        distance = np.abs(demandNodes[i] - servers).sum(axis=1)
        PreList[i] = np.argsort(distance)
    return PreList


def TransitionRate_divided(f, Pre_L, N, J, Lambda, Mu, n, nodeList):
    stateNum = 2 ** N
    units = list(range(N))
    up_row = np.zeros(n * len(nodeList))
    up_col = np.zeros(n * len(nodeList))
    up_rate = np.zeros(n * len(nodeList))
    down_row = np.zeros((N - n) * len(nodeList))
    down_col = np.zeros((N - n) * len(nodeList))
    down_rate = np.zeros((N - n) * len(nodeList))

    for i in range(len(nodeList)):
        busy = [m for m in range(N) if nodeList[i] & (1 << m) != 0]
        free = set(units) - set(busy)
        up_col[i * n: (i + 1) * n] = i
        up_row[i * n: (i + 1) * n] = nodeList[i] & ~(1 << np.array(busy))
        for j in range(J):
            for m in range(n):
                if Pre_L[j][m] in free:
                    break
                else:
                    up_rate[i * n + busy.index(Pre_L[j][m])] += f[j]

        down_row[i * (N - n): (i + 1) * (N - n)] = nodeList[i] | (1 << np.array(list(free)))
        down_col[i * (N - n): (i + 1) * (N - n)] = i
        down_rate[i * (N - n): (i + 1) * (N - n)] = Mu

    transup_rate = sparse.csc_matrix((up_rate * Lambda, (up_row, up_col)), shape=(stateNum, len(nodeList)))
    transdown_rate = sparse.csc_matrix((down_rate, (down_row, down_col)), shape=(stateNum, len(nodeList)))

    return transup_rate, transdown_rate


def BnD_divided(f, Pre_L, N, J, Lambda, Mu, p_n_B_new, n, nodeList):
    rho = Lambda / Mu
    transup, transdown = TransitionRate_divided(f, Pre_L, N, J, Lambda, Mu, n, nodeList)
    p_n_B_modify = (p_n_B_new * transup * n / rho + p_n_B_new * transdown * rho / (n + 1)) / (n * Mu + Lambda)

    del transup, transdown, p_n_B_new
    return [nodeList, p_n_B_modify]


def BnD_LargeScale(f, Pre_L, N, J, Lambda, Mu, batchsize, p_n_B_new, layer, processorsNum):
    # task assignment
    taskAssign = []
    minNum = batchsize
    print(minNum)
    for n in range(1, N):
        taskAssign.append([])
        taskNum = math.ceil(len(layer[n]) / minNum)
        if taskNum == 1:
            taskAssign[n - 1].append(layer[n])
        else:
            for j in range(taskNum):
                taskAssign[n - 1].append(layer[n][j:: taskNum])

    ctx = mp.get_context()
    ctx.reducer = Pickle4Reducer()
    pool = mp.Pool(processes=processorsNum)
    ite = 0
    time_list = []
    Err_list = []
    start_time = time.time()

    # start iteration
    while (ite < 1):
        p_n_B = np.copy(p_n_B_new)
        start_time_ite = time.time()
        for n in range(1, N):
            BnD_divided_partial = partial(BnD_divided, f, Pre_L, N, J, Lambda, Mu, p_n_B_new, n)
            new = pool.imap_unordered(BnD_divided_partial, taskAssign[n - 1])
            for i in new:
                p_n_B_new[i[0].astype(int)] = i[1]
        time_list += [time.time() - start_time_ite]
        ite += 1
        Err_list += [max(np.abs(p_n_B_new - p_n_B) / p_n_B_new)]
        print(ite)

    calcuTime = time.time() - start_time
    print("------ %s seconds, %s iterations ------" % (time.time() - start_time, ite))

    pool.close()
    pool.join()

    return np.array(time_list), np.array(Err_list), calcuTime


if __name__ == '__main__':
    seed = 4
    Data = {}
    Data['N'] = 26
    N = Data['N']
    Data['J'] = 71
    J = Data['J']
    Data['Mu'] = 1.76
    Mu = Data['Mu']
    Data['rho'] = 0.5
    rho = Data['rho']
    Data['Lambda'] = rho * Mu * N
    Lambda = Data['Lambda']
    Data['f'] = Random_Fraction(J, seed=seed)
    f = Data['f']
    Data['Pre_L'] = Random_Pref(N, J, seed=seed)
    Pre_L = Data['Pre_L']

    print('Number of units:', N)
    start_time = time.time()
    # steady state probability of busy units
    rho = Lambda / Mu
    p_n = np.array([(rho ** j) / math.factorial(j) for j in range(N + 1)])
    p_n = p_n / sum(p_n)

    Table = pd.DataFrame(columns=['processNum', 'time_list', 'Err_list', 'CalculationTime'])

    # initial states probabilities
    p_n_B = np.zeros(2 ** N)
    p_n_B_new = np.zeros((2, 2 ** N))
    for i in range(2 ** N):
        w = bin(i).count('1')
        p_n_B_new[:, i] = [w, 1 / comb(N, w)]

    layer = [np.sort(np.where(p_n_B_new[0, :] == i)[0]) for i in range(N + 1)]

    for i in range(2, 14, 2):
        processorsNum = i
        print(processorsNum)
        batchsize = 200000

        time_list, Err_list, calcuTime = BnD_LargeScale(f, Pre_L, N, J, Lambda, Mu, batchsize, p_n_B_new[1, :], layer,
                                                        processorsNum)

        Table.loc[len(Table.index)] = [i, time_list, Err_list, calcuTime]
        del p_n_B

    Table.to_csv('Table-Parallel-%s Units.csv' % N, index=None)