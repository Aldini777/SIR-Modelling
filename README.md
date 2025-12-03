#  SIR Modeling & Agent-Based Simulation Challenge

This repository contains a comprehensive analysis of epidemic dynamics, ranging from mathematical modeling using Differential Equations (**R**) to stochastic Agent-Based Simulations (**Python**).

The project is divided into 4 phases, exploring the behavior of infectious diseases under different scenarios such as vaccination, vital dynamics, and spatial distribution.

---

##  Important: How to View Graphs and Animations

To properly visualize all the interactive graphs and animations in this project, please follow these instructions:

1.  **Part 1 (R - Basic SIR):**
    * View `Reto_Parte_1.md` directly here on GitHub.
    * *(Optional)* Download `Reto_Parte_1.Rmd` and run it locally in RStudio.

2.  **Part 2 (R - Advanced Dynamics):**
    * View `Reto_Parte_2.md` for the static report and code explanation.
    * **For animations:** To view the time-lapse evolution (400 years), you must download `Reto_Parte_2.Rmd` and run it locally (Knit to HTML).

3.  **Parts 3 & 4 (Python - Spatial Simulation):**
    * These parts contain heavy animations using **Plotly**. GitHub cannot render them directly.
    * **Option A (Recommended):** Open the notebook directly in Google Colab to see the interactive bubbles:
        [** Open in Google Colab**](https://drive.google.com/file/d/1WC8bTOes5TTtwvk9c8tW8nDnLXNvmGFE/view?usp=sharing)
    * **Option B:** Download `Reto_Parte_3_y_4.ipynb` and run it locally using Jupyter Notebook or VS Code.

---

##  Project Structure

```text
SIR-Modelling/
├── Reto_Parte_1.md          # Static report for Phase 1
├── Reto_Parte_1.Rmd         # Source code (R) for Phase 1
├── Reto_Parte_2.md          # Static report for Phase 2
├── Reto_Parte_2.Rmd         # Source code (R) for Phase 2 (Includes Animations)
├── Reto_Parte_3_y_4.ipynb   # Python Notebook with Agent-Based Simulations
├── simulacion_ciudad.html   # Interactive HTML export of the simulation
└── README.md                # This file
```
**Phase 1: The Basic SIR Model (R)**
In this phase, we implemented the classic SIR (Susceptible-Infected-Recovered) model using Ordinary Differential Equations (ODEs).

Key Concepts:

$$R_0$$ (Basic Reproduction Number): Calculated as β/γ.

Threshold Theorem: Analytical derivation of the maximum number of infected individuals.

*S(t) Analysis*: Analytical approximation for the time required for Susceptibles to drop to half (N/2).

**Phase 2: Vital Dynamics & Vaccination (R)**
We expanded the model to include demographic factors and public health interventions.

1. Vital Dynamics (Births & Deaths)
We modified the differential equations to include birth (b) and death (μ) rates.

Short term (400 days): The epidemic behaves similarly to the standard SIR.

Long term (400 years): The model shows damped oscillations converging to an Endemic Equilibrium.

2. Vaccination Strategy
We analyzed the concept of Herd Immunity.
*Simulation*: We proved that vaccinating 80% of the population (when the threshold was 75%) successfully prevents an outbreak.

**Phase 3 & 4:** Agent-Based Simulation (Python)
Moving away from differential equations, we built a stochastic simulation where each individual is an autonomous agent moving in a 2D space.

Technologies: Python, NumPy, Pandas, Plotly Express.

**Scenarios Simulated:**
Square City:

Random movement with normal distribution.

Boundaries using "hard walls" logic.

Infection based on proximity radius (r).

Circular City:

Uniform distribution within a circle 

Movement constrained within radius D/2.

Cluster Distribution:

Simulates a "hotspot" (e.g., a concert or city center).

Agents are initialized using a Gaussian distribution around a random center.

Result: Extremely fast initial infection rates due to high local density.

Comparison of Variants:
We ran simulations comparing 4 scenarios for each geometry:

Base Case: Standard parameters.

High Density: Double population (N=600).

Social Distancing: Reduced infection radius (r=3.0).

Fast Recovery: Increased recovery rate (γ=0.2).

**Requirements**
To run the code locally, you need:

**For R:**
install.packages(c("deSolve", "ggplot2", "reshape2", "plotly"))
**For Python:**
pip install numpy pandas plotly

**Authors**
Aldo Resendiz Cravioto

Aldo Radamés Corral Verdugo

Virginia Díaz Moreno

November 2025
