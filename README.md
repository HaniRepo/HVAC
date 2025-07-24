
# HVAC System Simulator and MPC Evaluation (MATLAB)

This folder contains a **MATLAB-based simulator** and **evaluation scripts** built from real-world data to analyze and compare different system identification and predictive control models.  
It is part of the work described in the paper:  
**“A Switching Event-Triggered Predictive Control of HVAC Systems”** (also included in this repository).

---

##  Overview

-  A **simulator** is provided, built with real system data, to model HVAC dynamics.
-  An **automated script** evaluates different data segments (days, weeks, months) and compares:
  - Output of various **system identification models**.
  - Output of different **Model Predictive Controllers (MPC)**.
-  The purpose is to identify which MPC model captures **more realistic behavior** under different operating scenarios.

Because the original data is **confidential**, the repository includes **randomized demonstration data** to showcase the procedure.

---

##  Experimental Context

- Segments of experiments in this repository are referenced in the paper:
  > *A Switching Event-Triggered Predictive Control of HVAC Systems*
  
- A more detailed paper with expanded results is currently under review and is expected to be published by the **end of 2025**.

---

## 🗂 Repository Contents

- `setup.slx` – MATLAB Live Editor file with the simulation and evaluation setup.
- `setup.pdf` – PDF version of the Live Editor file for quick reference.
- `script.m` (or equivalent) – Automated script to:
  - Load dataset segments (day/week/month),
  - Run identification models,
  - Run MPC models,
  - Compare outputs against real (or demo) data.

---

##  How to Run

1. Ensure you have **MATLAB** installed with Simulink support.
2. Place the **data files** (provided demo data) in the same directory as the scripts.
3. Open `setup.slx` in MATLAB.
4. Simply run the script or Live Editor to execute the procedure.

>  **Tip:** You can also open the provided `setup.pdf` to get an overview before running the code.

---

##  Notes on Data

- Real data is **not included** due to confidentiality.
- A **random dataset** is provided to demonstrate the workflow and confirm that the system runs end-to-end.

---

##  References

If you use or build on this work, please cite:
> *A Switching Event-Triggered Predictive Control of HVAC Systems*  
(Paper available in the repository.)

---

## ✨ Future Work

A more elaborated version of the experimental setup and additional results will be published in a forthcoming paper by **end of 2025**.

---

## 🛠️ Requirements

- MATLAB (R2021b or newer recommended)
- Simulink support enabled

---

Enjoy experimenting and feel free to reach out for collaboration or questions!
