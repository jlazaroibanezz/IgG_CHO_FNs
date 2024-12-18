# IgG_CHO_FNs

This repository contains the code utilized to run the simulations performed in the paper "Multi-Scale Design and Optimization of Antibody Production via Flexible Nets" submitted to Computational and Structural Biotechnology Journal.

### 1. `CHOBiorFN_all_aa_nutrients_HP.py`

This script performs simulations to optimize antibody production in CHO cells. It models the cell culture conditions and tests different parameters to improve productivity. The simulation evaluates factors such as nutrient concentrations, nutrient uptake rates, dilution rates and growth rates to find the optimal conditions that maximize antibody yield.
The nutrient uptake rates were selected from Carinhas et al. 2013 and refer to the HP dataset mentioned in the paper.

### Key Features:

- Simulates various culture conditions and parameters.
- Uses optimization algorithms to improve antibody production.
- Outputs results for the maximum monoclonal antibody yield depending on nutrient concentrations, dilution rates, nutrient uptake rates and growth rates.

**How to use**: Simply run the script, and it will automatically simulate and return the optimized antibody flux depending on the previous parameters, among other interesting values.

```
$ python3 CHOBiorFN_all_aa_nutrients_HP.py
```

### 2. `nutrients_HP.py`

This script contains four main functions:

1. Limiting Metabolite Discovery (genCHOBiorFNmax and genCHOBiorFNmin): Identify the metabolite that limits the growth and antibody production in the cell culture.
2. Culture Medium Minimization (MinimizeMediumim and EcoMinimizeMediumim): MinimizeMediumim implements a minimization algorithm to reduce the number of components in the culture medium while maintaining antibody productivity while EcoMinimizeMediumim minimizes an economically ponderated sum of the nutrient concentrations based on the price of each nutrient in the market.

### Key Features:

- Analyzes cell culture data to identify limiting metabolites.
- Minimizes the culture medium formulation without sacrificing productivity.
- Outputs the optimal medium composition and the limiting metabolites that need attention for enhanced productivity.

**How to use**: This script is automatically called by CHOBiorFN_all_aa_nutrients_HP.py 

### Requirements

Both scripts require the following Python libraries:

- numpy
- cobra
- pandas
- fnyzer
- pandas_ods_reader

and the previous installation of a a solver that will solve the linear programming problem. We suggest [GLPK](https://www.gnu.org/software/glpk/), [Gurobi](https://www.gurobi.com/) or
CPLEX. In our case, we used [CPLEX](https://www.ibm.com/es-es/products/ilog-cplex-optimization-studio).
