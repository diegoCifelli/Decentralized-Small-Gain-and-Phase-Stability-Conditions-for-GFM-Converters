# Code for “Decentralized Small Gain and Phase Stability Conditions for Grid-Forming Converters: Limitations and Extensions”

This repository contains the scripts, data, and supplementary materials used in the paper:

> **“Decentralized Small Gain and Phase Stability Conditions for Grid-Forming Converters: Limitations and Extensions”**  
> Authors: Diego Cifelli, Adolfo Anta  
> Submitted to PCSS 2026.

---

## 📘 Overview

This repository provides all the necessary files to reproduce the numerical results and figures presented in the paper.  
The scripts are organized to allow straightforward replication of the study and exploration of alternative scenarios.

---

## ⚙️ How to Use

The scripts in this repository are based on the **Simplus Tool**: [https://github.com/Future-Power-Networks/Simplus-Grid-Tool](https://github.com/Future-Power-Networks/Simplus-Grid-Tool).  
Before using any of the scripts in this repository, please download and install the Simplus Tool.

The GFM converter model used in the paper is provided in the file `GridFormingVSI.m`, which should be copied and used to **overwrite** the original file located in the folder `+SimplusGT/+Class`.

Three models are provided as Excel sheets:
- `UserData_inf_bus.xlsm`: Model of the GFM converter connected to an infinite bus.  
- `UserData_inf_bus_Fig_8.xlsm`: Model of the GFM converter connected to an infinite bus with reduced damping, used to reproduce the unstable case in *Example 1* of the paper.  
- `UserData_IEEE14_Unstable.xlsm`: Model of the IEEE 14-bus system, used in *Example 2* of the paper.

The main script is `UserMain.m`, which reads the model defined in one of the Excel files and generates all the necessary objects for the analysis. This script must be run first.

After that:
- For the infinite bus cases, run `dec_conditions_inf_bus.m`.  
- For the IEEE 14-bus case, run `dec_conditions_multibus.m`.

## 📫 Contact

For questions or further information, please contact:  
**Diego Cifelli**  
📧 [diego.cifelli@ait.ac.at]
