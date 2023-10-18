# -*- coding: utf-8 -*-

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

# Custom ForkingPickler for multiprocessing with protocol 4
class ForkingPickler4(ForkingPickler):
    # Custom Pickler class derived from ForkingPickler
    def __init__(self, *args):
        if len(args) > 1:
            args[1] = 2
        else:
            args.append(2)
        super().__init__(*args)

    @classmethod
    def dumps(cls, obj, protocol=4):
        # Custom method to pickle an object with the specified protocol
        return ForkingPickler.dumps(obj, protocol)

def dump(obj, file, protocol=4):
    # Function to pickle an object and write it to a file
    ForkingPickler4(file, protocol).dump(obj)

class Pickle4Reducer(AbstractReducer):
    # Custom Pickle Reducer class derived from AbstractReducer
    ForkingPickler = ForkingPickler4
    register = ForkingPickler4.register
    dump = dump

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

    # Initialize arrays to store preference lists and distances of atom-server pairs
    PreList = np.zeros((J, N))

    # Loop over geographical atoms
    for i in range(J):

        # Calculate the Manhattan distances between the current geographical atom and all servers
        distance = np.abs(demandNodes[i] - servers).sum(axis=1)

        # Generate the preference list of the servers based on the distances from near to far
        PreList[i] = np.argsort(distance)
        
    return PreList

# Function to generate transition rates for a subset of states
def TransitionRate_divided(f, Pre_L, N, J, Lambda, Mu, n, nodeList):
    '''
    :param f: fraction of region-wide demand
    :param Pre_L: preference list
    :param N: Number of units
    :param J: Number of geographical atoms
    :param Lambda: overall arrival rate within the entire region
    :param Mu: service rate of each response unit 
    :param n: current layer
    :param nodeList: states need to generate the transition rates
    :return: transition rate matrix of the divided system with a subset of states
    '''
    stateNum = 2 ** N
    units = list(range(N))

    # Initialize arrays to store the upward and downward transition rates
    up_row = np.zeros(n * len(nodeList))
    up_col = np.zeros(n * len(nodeList))
    up_rate = np.zeros(n * len(nodeList))
    down_row = np.zeros((N - n) * len(nodeList))
    down_col = np.zeros((N - n) * len(nodeList))
    down_rate = np.zeros((N - n) * len(nodeList))

    # Loop over all states in the nodesList
    for i in range(len(nodeList)):
        # Identify busy and free units in the current state
        busy = [m for m in range(N) if nodeList[i] & (1 << m) != 0]
        free = set(units) - set(busy)

        # Populate column i of the upward transition matrix
        up_col[i * n: (i + 1) * n] = i

        # Find all unit-step reachable states by changing exactly one unit from 1 to 0 of the current state
        up_row[i * n: (i + 1) * n] = nodeList[i] & ~(1 << np.array(busy))

        # Loop over geographical atoms
        for j in range(J):

            # For each atom, find the optimal dispatch unit and add the fraction of region-wide demand to the corresponding location in upward transition matrix
            for m in range(n):
                if Pre_L[j][m] in free:
                    break
                else:
                    up_rate[i * n + busy.index(Pre_L[j][m])] += f[j]

        # Populate column i of the downward transition matrix
        # Find all unit-step reachable states by changing exactly one unit from 0 to 1 of the current state
        down_row[i * (N - n): (i + 1) * (N - n)] = nodeList[i] | (1 << np.array(list(free)))
        down_col[i * (N - n): (i + 1) * (N - n)] = i

        # Fill ``Mu''  to corresponding locations in the downward transition matrix
        down_rate[i * (N - n): (i + 1) * (N - n)] = Mu

    # Create sparse matrices for upward and downward transition rates
    transup_rate = sparse.coo_matrix((up_rate * Lambda, (up_row, up_col)), shape=(stateNum, len(nodeList)))
    transdown_rate = sparse.coo_matrix((down_rate, (down_row, down_col)), shape=(stateNum, len(nodeList)))

    return transup_rate, transdown_rate

# Function to solve the birth and death model for a subset of states
def BnD_divided(f, Pre_L, N, J, Lambda, Mu, p_n_B_new, n, nodeList):
    '''
    :param f: fraction of region-wide demand
    :param Pre_L: preference list
    :param N: Number of units
    :param J: Number of geographical atoms
    :param Lambda: overall arrival rate within entire region
    :param Mu: service rate of each response unit 
    :param p_n_B_new: last computed conditional probability of the subset of states
    :param n: current layer
    :param nodeList: states need to generate the transition rates
    :return: node list and its updated probability distribution solved by BnD model 
    '''
    rho = Lambda / Mu

    # Calculate transition rates for the divided system consisting of the subset of states
    transup, transdown = TransitionRate_divided(f, Pre_L, N, J, Lambda, Mu, n, nodeList)

    # Update the probabilities using the birth and death model equations
    p_n_B_modify = (p_n_B_new * transup * n / rho + p_n_B_new * transdown * rho / (n + 1)) / (n * Mu + Lambda)

    del transup, transdown, p_n_B_new
    return [nodeList, p_n_B_modify]

# Function to perform BnD iteration for a large-scale problem using parallel computing
def BnD_LargeScale(f, Pre_L, N, J, Lambda, Mu, batchsize, p_n_B_new, layer, processorsNum):
    '''
    :param f: fraction of region-wide demand
    :param Pre_L: preference list
    :param N: Number of units
    :param J: Number of geographical atoms
    :param Lambda: overall arrival rate within the entire region
    :param Mu: service rate of each response unit 
    :param batchsize: number of states in each task
    :param p_n_B_new: initial conditional probability distribution
    :param layer: the set of states with a certain number of busy units
    :param processorsNum: total number of processors
    :return: probability distribution solved by BnD model and its execution time  
    '''
    
    # Task assignment
    taskAssign = []
    minNum = batchsize
    print(minNum)
    for n in range(1, N):
        taskAssign.append([])

        # Total number of batched tasks
        taskNum = math.ceil(len(layer[n]) / minNum)
        if taskNum == 1:
            taskAssign[n - 1].append(layer[n])
        else:
            for j in range(taskNum):
                # Assign s(=batchsize) states to each tasks
                taskAssign[n - 1].append(layer[n][j:: taskNum])

    # Use protocol 4 (Not necessary when batchsize is small)
    ctx = mp.get_context()
    ctx.reducer = Pickle4Reducer()

    # Create a pool with a given number of processors
    pool = mp.Pool(processes=processorsNum)

    # Initialize variables for iteration
    ite = 0
    time_list = []
    Err_list = []
    start_time = time.time()

    # Start iteration
    while (ite < 1):
    # Here we only run a single iteration. To run the full iterations, unmark the following command:
    # while gap >= 1e-3:
        
        p_n_B = np.copy(p_n_B_new)
        start_time_ite = time.time()

        # Iterate over layers using updating equations in the Birth and Death model
        for n in range(1, N):
            # Define the partial function with all parameters except the node list fixed to perform the updating for tasks
            BnD_divided_partial = partial(BnD_divided, f, Pre_L, N, J, Lambda, Mu, p_n_B_new, n)

            # Divide the tasks in the current layer to processors with task assignment list
            new = pool.imap_unordered(BnD_divided_partial, taskAssign[n - 1])

            # Conquer the results from processors 
            for i in new:
                p_n_B_new[i[0].astype(int)] = i[1]

        # Update iteration count and calculate the maximum relative change
        time_list += [time.time() - start_time_ite]
        ite += 1
        Err_list += [max(np.abs(p_n_B_new - p_n_B) / p_n_B_new)]
        print(ite)

    # Record the total calculation time
    calcuTime = time.time() - start_time
    print("------ %s seconds, %s iterations ------" % (time.time() - start_time, ite))

    # Close the processor pool
    pool.close()
    pool.join()

    return p_n_B_new, np.array(time_list), np.array(Err_list), calcuTime


if __name__ == '__main__':
    # Set a random seed for reproducibility
    seed = 4

    # Initialize a dictionary to store critical data
    Data = {}

    # Set the number of server units
    Data['N'] = 26
    N = Data['N']

    # Set the number of geographical atoms
    Data['J'] = 71
    J = Data['J']

    # Set the service rate
    Data['Mu'] = 1.76
    Mu = Data['Mu']

    # # Set the utilization and overall arrival rate within the entire region
    Data['rho'] = 0.5
    rho = Data['rho']
    Data['Lambda'] = rho * Mu * N
    Lambda = Data['Lambda']

    # Generate a random fraction of region-wide demand
    Data['f'] = Random_Fraction(J, seed=seed)
    f = Data['f']

    # Generate random preference list and distance list
    Data['Pre_L'] = Random_Pref(N, J, seed=seed)
    Pre_L = Data['Pre_L']

    print('Number of units:', N)
    start_time = time.time()
    
    # Record results in the DataFrame
    Table = pd.DataFrame(columns=['processNum', 'time_list', 'Err_list', 'CalculationTime'])

    # Calculate the steady state probability of busy units
    rho = Lambda / Mu
    p_n = np.array([(rho ** j) / math.factorial(j) for j in range(N + 1)])
    p_n = p_n / sum(p_n)
    
    # Initialize states probabilities by assuming all states within a layer is equally distributed
    p_n_B_new = np.zeros((2, 2 ** N))
    for i in range(2 ** N):
        w = bin(i).count('1')
        p_n_B_new[:, i] = [w, 1 / comb(N, w)]

    # Define layers based on the number of busy units
    layer = [np.sort(np.where(p_n_B_new[0, :] == i)[0]) for i in range(N + 1)]

    # Loop over number of processors
    for i in range(2, 14, 2):
        # Set the number of processors
        processorsNum = i
        print(processorsNum)

        # Set the number of batchsize
        batchsize = 200000

        # Calculate the steady state conditional probability with parallel computing
        p_n_B_new, time_list, Err_list, calcuTime = BnD_LargeScale(f, Pre_L, N, J, Lambda, Mu, batchsize, p_n_B_new[1, :], layer, processorsNum)

        # Recover the probability distribution
        for i in range(N+1):
            p_n_B_new[1, layer[i]] =  p_n[i] * p_n_B_new[1, layer[i]]

        # Record results in the DataFrame
        Table.loc[len(Table.index)] = [i, time_list, Err_list, calcuTime]

        # Clean up variables and release the memory for next iteration
        del p_n_B_new

    # Save the results to a CSV file
    Table.to_csv('Table-Parallel-%s Units.csv' % N, index=None)
