# IgG_CHO_FNs

This repository contains the code utilized to run the simulations performed in the paper "Flexible Nets to optimize antibody production in Chinese Hamster Ovary cells" submitted to Computational and Structural Biotechnology Journal.

## 1. CHOBiorFN_all_aa_nutrients_HP_git.py

This script performs simulations to optimize antibody production in CHO cells. It models the cell culture conditions and tests different parameters to improve productivity. The simulation evaluates factors such as nutrient concentrations, nutrient uptake rates, dilution rates and growth rates to find the optimal conditions that maximize antibody yield.
The nutrient uptake rates were selected from Carinhas et al. 2013 and refer to the HP dataset mentioned in the paper.

### Key Features:

- Simulates various culture conditions and parameters.
- Uses optimization algorithms to improve antibody production.
- Outputs results for the best maximum monoclonal antibody yield depending on nutrient concentrations, dilution rates, nutrient uptake rates and growth rates.
