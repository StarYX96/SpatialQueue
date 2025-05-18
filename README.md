# Replication Package for *A Geometrically Convergent Solution to Spatial Hypercube Queueing Models*
**By Cheng Hua, Jun Luo, Arthur J. Swersey, Yixing Wen**  

In the main context of our paper, we evaluate the performance of our proposed Conditional Probability Update (CPU) algorithms as well as its parallel version. This ReadMe file provides detailed instructions for replicating the figures, experiments, and analyses presented in the paper. For any issues or questions, please contact jluo_ms@sjtu.edu.cn.

## 1. General Introduction  
The replication package includes Python code files to reproduce all figures and tables in the paper. Requirements:  
- Python >3.9  
- Install dependencies:  
```bash
pip install -r requirements.txt
```

## 2. Experiment Settings and Code Files  

### Table 1: List of Experiment Code Files  
| Results                  | File Names of Codes                          | Introduction                                                                 |
|--------------------------|---------------------------------------------|-----------------------------------------------------------------------------|
| CPU - Heterogeneous Rate | `HeterogeneousRate-Greenville.py`<br>`HeterogeneousRate-StPaul.py`<br>`Figure 5.py` | Experiments under heterogeneous rate cases [1] |
| CPU - Homogeneous Rate   | `HomogeneousRate-Greenville.py`<br>`HomogeneousRate-StPaul.py`<br>`Figure 6.py`   | Experiments under homogeneous rate cases [1]   |
| Discrete Event Simulation| `Table3-Simulation.py`                      | St.Paul simulation with different service time distributions [1]            |
| Parallel CPU             | `Table4-Parallel.py`                        | Parallel computing experiments [1]            |

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
