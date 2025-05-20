# Replication Package for *A Geometrically Convergent Solution to Spatial Hypercube Queueing Models*
**By Cheng Hua, Jun Luo, Arthur J. Swersey, Yixing Wen**  

In the main context of our paper, we evaluate the performance of our proposed Conditional Probability Update (CPU) algorithms as well as its parallel version. This file provides detailed instructions for replicating the figures, experiments, and analyses presented in the paper. For any issues or questions, please contact jluo_ms@sjtu.edu.cn.

## 1. General Introduction  
The replication package includes Python code files, structured to facilitate the reproduction of all figures and tables in the paper. The package is compatible with both Linux and Windows systems. To compile and run the codes successfully, please use Python with a version higher than 3.9. Also, to produce the results, users need to prepare their Python environment by installing the following file:  
```bash
pip install -r requirements.txt
```
In what follows, we outline the experiment configurations and provide a comprehensive guide to the code files.

## 2. Experiment Settings and Code Files  
The paper tests various scenarios by varying the number of units and system utilization. Numerical experiments are conducted to report the computation time, mean response time (MRT), unit utilization, and the maximum percentage relative error (MPRE) of different measures between methods. The test experiment hardware setup used by the author is as follows:
-	**CPU**: AMD Ryzen 9 5950X 16-Core Processor with Nvidia GeForce RTX 3090 Graphics
-	**Memory (RAM)**: 64GB
-	**Storage**: 1.1TB

Notice that the computation time of experiments are corelated to machines’ execution status (e.g., temperature, available memory, version of operation system, etc.). So it’s common that the computation time might be slightly different.

Except for the discrete event simulation procedure, experiments do not need to replicate since they take the same distribution of units’ location and fraction of calls under the same scenarios. The code files of experiments of two cities lie solely in the city’s structure and system’s estimated parameters, with no changes in algorithm design. All experiments will record their computation time and MPRE of certain measures comparing to other approaches. Please refer to our paper for further details.  The experiments are listed in Table 1.

For each experiment, the result is exported to a CSV file with the corresponding filename after execution. It is important to note that the tables presented in the paper have undergone formatting adjustments to enhance readability, such as unifying scientific notation scales and highlighting specific results in bold. These formatting adjustments require post-processing. The tables can be directly reproduced by the above code files representatively. 


## Table 1 List of Experiment Code Files and Corresponding Introductions
<table style="width: 100%; table-layout: fixed;margin: 0 auto;">
    <colgroup>
    <col style="width: 30%;">
    <col style="width: 30%;">
    <col style="width: 40%;">
</colgroup>
    <tr>
        <th style="text-align: center"> Results </th>
        <th style="text-align: center"> File Names of Codes </th>
        <th style="text-align: center; vertical-align: middle"> Introduction </th>
    </tr>
    <tr>
        <td colspan="3" align="center">CPU - Heterogeneous Rate</td>
    </tr>
    <tr>
      <td rowspan="3">Figure 5</td>
      <td>HeterogeneousRate- StPaul.py</td>
      <td rowspan="3" style="width: 40%;">The experiments of two cities under heterogeneous rate cases.  The third file reproduces Figure 5, which relies on the CSV results of corresponding experiments given by the first two files. For brevity, the CSV results are provided in prior. </td>
    </tr>
    <tr>
      <td>HeterogeneousRate-Greenville.py</td>
    </tr>
    <tr>
      <td>Figure 5.py</td>
    </tr>
     <tr>
        <td colspan="3" align="center">CPU - Homogeneous Rate Cases</td>
    </tr>
    <tr>
      <td rowspan="3">Figure 6</td>
      <td>HomogeneousRate- StPaul.py</td>
      <td rowspan="3" style="text-align: justify; vertical-align: top; word-wrap: break-word">The experiments of two cities under homogeneous rate cases.  The third file reproduces Figure 6, which relies on the CSV results of corresponding experiments given by the first two files. For brevity, the CSV results are provided in prior. </td>
    </tr>
    <tr>
      <td>HomogeneousRate-Greenville.py</td>
    </tr>
    <tr>
      <td>Figure 6.py</td>
    </tr>
    <tr>
        <td colspan="3" align="center">Discrete Event Simulation</td>
    </tr>
    <tr>
        <td>Table 3</td>
        <td>Table3-Simulation.py</td>
        <td style="text-align: justify; vertical-align: top; word-wrap: break-word">The experiments of discrete event simulation of St. Paul. under heterogeneous rate cases and different service time distribution.</td>
    </tr>
    <tr>
        <td colspan="3" align="center">Parallel CPU</td>
    </tr>
    <tr>
        <td>Table 4</td>
        <td>Table4-Parallel.py</td>
        <td style="text-align: justify; vertical-align: top; word-wrap: break-word">The experiments of parallel experiment of St. Paul. under homogeneous rate cases.</td>
    </tr>
</table>
