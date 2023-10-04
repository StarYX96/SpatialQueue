# BirthDeathChain-SpatialQueue

On this site, we provide the codes used to conduct the numerical experiments of the paper entitled "A Birth and Death Chain Solution to a Spatial Queueing Problem". All the cosed are written in Python. To compile and run the codes successfully, one may need to install Python with veirsion larger than 3.9.

---File HomogeneousRate.py is the Python codes that calculates the probability distribution for homogeneoues service rate using both original hypercube model and 
the birth and death chain model. This code will output the coefficient generation time, probabilities updating time and maximum percentage relative error (MPRE) of both methods, which can be used to calculate the Table 2 and Table EC.4.3. See Section 7.1 of the paper for detailed post-processing method.

---File HeterogeneousRate.py is the Python codes that calculates the probability distribution for heterogeneoues service rate using both sparse solver and 
the birth and death chain model. This code will output the probability updating time and MPRE of both methods.

---File Parallel.py is the Python codes that calculates the probability distribution using parallel computing framework. This code will output

---Directory "Codes for EC" contains all the Python codes that are used to conduct the numerical experiments listed in Section EC.4 of the paper entitled "A Birth and Death Chain Solution to a Spatial Queueing Problem".
```diff
- Files and descriptions of them.
```
 

