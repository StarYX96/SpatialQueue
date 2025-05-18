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

### <center>Table 1 List of Experiment Code Files and Corresponding Introductions</center>
| Results  | File Names of Codes | Introduction |
|--------------------------|---------------------------------------------|-----------------------------------------------------------------------------|
| CPU - Heterogeneous Rate | `HeterogeneousRate-Greenville.py`<br>`HeterogeneousRate-StPaul.py`<br>`Figure 5.py` | Experiments under heterogeneous rate cases [1] |
| CPU - Homogeneous Rate   | `HomogeneousRate-Greenville.py`<br>`HomogeneousRate-StPaul.py`<br>`Figure 6.py`   | Experiments under homogeneous rate cases [1]   |
| Discrete Event Simulation| `Table3-Simulation.py`                      | St.Paul simulation with different service time distributions [1]            |
| Parallel CPU             | `Table4-Parallel.py`                        | Parallel computing experiments [1]            |

<table>
  <thead>
    <tr>
      <th rowspan="2">Results</th>
      <th colspan="2">File Names of Codes</th>
      <th rowspan="2">Introduction</th>
    </tr>
    <tr>
      <th>City Implementation</th>
      <th>Visualization</th>
    </tr>
  </thead>
  <tbody>
    <!-- 异构速率场景 -->
    <tr>
      <td rowspan="2">CPU - Heterogeneous Rate</td>
      <td>HeterogeneousRate-Greenville.py</td>
      <td rowspan="2">Figure 5.py</td>
      <td rowspan="2">Two cities experiments under heterogeneous rates [1]</td>
    </tr>
    <tr>
      <td>HeterogeneousRate-StPaul.py</td>
    </tr>

    <!-- 同构速率场景 -->
    <tr>
      <td rowspan="2">CPU - Homogeneous Rate</td>
      <td>HomogeneousRate-Greenville.py</td>
      <td rowspan="2">Figure 6.py</td>
      <td rowspan="2">Two cities experiments under homogeneous rates [1]</td>
    </tr>
    <tr>
      <td>HomogeneousRate-StPaul.py</td>
    </tr>

    <!-- 独立实验行 -->
    <tr>
      <td>Discrete Event Simulation</td>
      <td colspan="2">Table3-Simulation.py</td>
      <td>St.Paul simulation with service time variations [1]</td>
    </tr>
    
    <tr>
      <td>Parallel CPU</td>
      <td colspan="2">Table4-Parallel.py</td>
      <td>Parallel computing experiments [1]</td>
    </tr>
  </tbody>
</table>

### File Naming Conventions:  
- City-specific implementations:  
  - `-Greenville.py`: Greenville city structure  
  - `-StPaul.py`: St.Paul city structure  
- Visualization scripts require CSV outputs from corresponding experiments  

## 3. Hardware Specifications  
```markdown
- CPU: AMD Ryzen 9 5950X 16-Core Processor + NVIDIA RTX 3090
- RAM: 64GB DDR4
- Storage: 1.1TB NVMe SSD
```

## 4. Results Reproduction Notes  
> **Important**: Tables in the paper include post-processing adjustments:  
> - Unified scientific notation scales  
> -Bold highlighting of key results  
> Use scripts in `postprocessing/` directory to format raw CSV outputs [1].

## Contact Information  
For technical support contact:  
`jluo_ms@sjtu.edu.cn`  
