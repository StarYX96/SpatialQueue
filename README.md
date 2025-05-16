# Spatial Queue (Need Modified)

On this site, we provide the codes used to conduct the numerical experiments for the paper titled "A Geometrically Convergent Solution to Spatial Hypercube Queueing Models for Emergency Services". All codes are written in Python. To compile and run the codes successfully, please use Python with a version higher than 3.9.

---File HomogeneousRate.py is the Python code that calculates the probability distribution for homogeneous service rates using both the original hypercube method and our approach. This code outputs the coefficient generation time, probability updating time, and maximum percentage relative error (MPRE) of both methods, which can be used to generate Table 2 and Table EC.3.1. 

---File HeterogeneousRate.py is the Python code that calculates the probability distribution for heterogeneous service rates using both the sparse solver and our approach. This code outputs the probability updating time and MPRE of both methods, which can be used to generate Table 4 and Table EC.3.3.

---File Parallel.py is the Python code that calculates the probability distribution using a parallel computing framework. This code outputs the total computation time for different numbers of workers, which can be used to generate Table 6. Note that one may need to do the linear regression additionaly to calculate the $P$ and $R^2$.

---File requirements.txt lists all the modules (package and libraries) we use to execute the experiments.

<!--Directory "Codes for EC" contains all the Python codes that are used to conduct the numerical experiments listed in Section EC.3 of the paper entitled "A Geometrically Convergent Solution to Spatial Hypercube Queueing Models for Emergency Services".
   1) File HomogeneousRate-ModifiedHypercube.py is the Python code that calculates the probability distribution for homogeneoues service rate using modified original hypercube model. This code will output the coefficient generation time and probabilities updating time, which can be used to calculate the 'MOH' parts of Table EC.4.3 and Table EC.4.4. See Section EC.4.2 of the paper for detailed post-processing method.
   2) Files HomogeneousRate-Greenville.py, HeterogeneousRate-Greenville.py and HomogeneousRate-ModifiedHypercube.py are Python codes that calculates the probability distribution using the dataset of Greenville County, South Carolina. The contents of them are similar to what we introduced above. These codes help us to calculate TableS EC.4.1-EC.4.2, and Tables EC.4.6-EC.4.8-->

