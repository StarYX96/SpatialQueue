# BirthDeathChain-SpatialQueue

On this site, we provide the codes used to conduct the numerical experiments for the paper titled "A Birth and Death Chain Solution to a Spatial Queueing Problem." All codes are written in Python. To compile and run the codes successfully, one needs to install Python with a version higher than 3.9.

---File HomogeneousRate.py is the Python code that calculates the probability distribution for homogeneous service rates using both the original hypercube model and the birth and death chain model. This code outputs the coefficient generation time, probability updating time, and maximum percentage relative error (MPRE) of both methods, which can be used to generate Table 2, Table EC.4.3, and Table EC.4.4. 

---File HeterogeneousRate.py is the Python code that calculates the probability distribution for heterogeneous service rates using both the sparse solver and the birth and death chain model. This code outputs the probability updating time and MPRE of both methods, which can be used to generate Table 3 and Table EC.4.5.

---File Parallel.py is the Python code that calculates the probability distribution using a parallel computing framework. This code outputs the total computation time for different numbers of workers, which can be used to generate Table 4. Note that one may need to do the linear regression additionaly to calculate the $P$ and $R^2$.

---Directory "Codes for EC" contains all the Python codes that are used to conduct the numerical experiments listed in Section EC.4 of the paper entitled "A Birth and Death Chain Solution to a Spatial Queueing Problem".
   1) File HomogeneousRate-ModifiedHypercube.py is the Python code that calculates the probability distribution for homogeneoues service rate using modified original hypercube model. This code will output the coefficient generation time and probabilities updating time, which can be used to calculate the 'MOH' parts of Table EC.4.3 and Table EC.4.4. See Section EC.4.2 of the paper for detailed post-processing method.
   2) Files HomogeneousRate-Greenville.py, HeterogeneousRate-Greenville.py and HomogeneousRate-ModifiedHypercube.py are Python codes that calculates the probability distribution using the dataset of Greenville County, South Carolina. The contents of them are similar to what we introduced above. These codes help us to calculate TableS EC.4.1-EC.4.2, and Tables EC.4.6-EC.4.8
 

