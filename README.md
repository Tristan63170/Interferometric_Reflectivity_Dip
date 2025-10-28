# Interferometric Reflectivity Dip study

_This repository is dedicated to the study of "the Interferometric Reflectivity Dip" a new plasmonic effect predicted theoretically by the PHOTON/ELENA team from Institut Pascal (Clermont-Ferrand in France)._

## 📘 Overview

Surface Plasmon Resonance (SPR) sensors are used to measure biomolecular interactions in biological analyses. Current SPR sensors, which rely on measuring the reflection coefficient, are considered to be operating near the fundamental limit of sensitivity.

In 2016, a theoretical scientific paper suggested that SPR could induce a new effect—beam narrowing—which has the potential to push the accepted boundaries of SPR sensor sensitivity. This beam narrowing is thought to result from interference between the reflected beam (originating from the incoming beam's reflection on a multilayered structure) and the surface plasmon (excited within the structure by the incoming beam). This effect is now referred to as the Interferometric Reflectivity Dip.

The goal of this repository is to carry out a complete numerical analysis of the effect predicted in 2016, with the aim of proposing a realistic setup to observe it and evaluating its potential impact on SPR sensor sensitivity. The first step is to identify the experimental parameters that maximize the effect, then design a corresponding realistic setup, and finally compare its sensitivity to that of current SPR sensors.

Two versions of the study have been conducted, depending on the material at the end of the multilayered structure. In the first version, air is considered, corresponding to the 2016 paper. In the second version, water is used, which better reflects real-world SPR sensor conditions. Each version is stored in a dedicated folder; the code is largely similar across both, but the results differ.

The entire study has been conducted using Python. In particular, beam–multilayer structure interactions are simulated using the PyMoosh library. Most of the analysis is performed within Jupyter Notebooks, which include comments explaining much of the work carried out. Additional details about this repository are provided below.

This study is part of an ongoing research effort, and a related scientific paper is currently being written.

## 🗂 Repository Structure

```
.
├── my_module.py                      # Script with functions for beam and optic systems simulations
├── with_air/                         # Studies and analysis with air
│   ├── full_system.ipynb             # Notebook dedicated to the experimental setup determination
│   ├── observable_determination_and_sensitivity.ipynb  # Notebook dedicated to the sensitivity study
│   ├── optimal_parameters_determination.ipynb          # Notebook dedicated to the optimal parameters determination
│   ├── observable_determination_and_sensitivity_parallel.ipynb  # Notebook dedicated to the sensitivity study with multi-CPU computing (faster)
│   ├── optimal_parameters_determination_parallel.ipynb          # Notebook dedicated to the optimal parameters determination with multi CPU-computing (faster)
│   └── Plots/
│       ├── pdf/                      # PDF plots
│       └── svg/                      # SVG plots
├── with_water/                       # Studies and analysis with water
│   ├── animation_generator.py        # Script to generate animations of the effect
│   ├── full_system.ipynb             # Notebook dedicated to the experimental setup determination
│   ├── r_expression_study.ipynb  # Notebook dedicated to the theoretical study of the reflection coefficient
│   ├── surface_plasmon_illustration.ipynb  # Notebook dedicated to the production of 2D surface plasmon illustration
│   ├── observable_determination_and_sensitivity.ipynb  # Notebook dedicated to the sensitivity study
│   ├── optimal_parameters_determination.ipynb          # Notebook dedicated to the optimal parameters determination
│   ├── observable_determination_and_sensitivity_parallel.ipynb  # Notebook dedicated to the sensitivity study with multi-CPU computing (faster)
│   ├── optimal_parameters_determination_parallel.ipynb          # Notebook dedicated to the optimal parameters determination with multi CPU-computing (faster)
│   ├── NL_check.ipynb  # Notebook dedicated to the study of the non locality impact and to the production of pedagogical plots about observables
│   ├── H_module.mp4                  # Animation of the field deformation in log scale
│   ├── sqrt_H_module.mp4             # Animation of the square root of the field deformation
│   └── Plots/
│       ├── for_the_animation/        # Specific plots for animation (each frame)
│       ├── pdf/                      # PDF plots
│       └── svg/                      # SVG plots
```

## 🛠 Requirements

- Python >= 3.13 

Main dependencies:

- PyMoosh
- numpy  
- matplotlib  
- scipy  
- jupyter  
- imageio
- pickle
- joblib
- tqdm
- multiprocessing

## 👤 Authors

- Tristan Miralles (tristan.miralles@uca.fr)  
- Antoine Moreau
- Pauline Bennet
- Denis Langevin
