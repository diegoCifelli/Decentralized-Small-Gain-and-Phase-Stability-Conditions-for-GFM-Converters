# Code for “Paper Title”

This repository contains the scripts, data, and supplementary materials used in the paper:

> **“Paper Title”**  
> *Authors: Diego Cifelli, Adolfo Anta
> Submitted to PCSS 2026.

---

## 📘 Overview

This repository provides all the necessary files to reproduce the numerical results and figures presented in the paper.  
The scripts are organized to allow straightforward replication of the study and exploration of alternative scenarios.

---

## ⚙️ How to Use

The scripts are based on the Simplus Tool: https://github.com/Future-Power-Networks/Simplus-Grid-Tool.  
Before using the script in this reposory, download and install the Simplus Tool.

The GFM converter model used in the paper is provided the file GridFormingVSI.m, which need to be copy and overwrite the 
original file in the folder +SimplusGT/+Class.

Three models are provided as excel sheets:
- UserData_inf_bus.xlms: The model of the GFM converter connected to an infiinite bus.
- UserData_inf_bus_Fig_8.xlms: The model of the GFM converter connected to an infiinite bus, which denuted damping, to reproduce the unstable case in Example 1 of the paper.
- UserData_IEEE14_Unstable.xlms: The model of the IEEE14 Bus system, used in example 2 of the paper.

The script UserMain.m is the script that reads the model defined in one of the excel file, and generate all the object needed for the analysis. This is the first script that must be run.

For the infinite bus cases, then run the script "dec_conditions_inf_bus.m", for the IEEE14 case run "dec_conditions_multibus.m".


## 📫 Contact

For questions or further information, please contact:  
**Diego Cifelli**  
📧 [diego.cifelli@ait.ac.at]
