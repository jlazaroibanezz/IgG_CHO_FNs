# Genome-Scale Metabolic Model Simulations for CHO Cells

This repository provides Python scripts to work with the genome-scale metabolic model `iCHOv1`, simulating antibody production in Chinese Hamster Ovary (CHO) cells. The [iCHOv1](http://bigg.ucsd.edu/models/iCHOv1) model is a genome-scale metabolic model of Chinese hamster ovary (CHO) cells, widely used for therapeutic protein production. It includes 1766 genes, 6663 reactions, and 4456 metabolites, covering key pathways like glycolysis, TCA cycle, and amino acid biosynthesis. The model supports optimizing antibody production, analyzing metabolic fluxes, and improving bioprocess conditions.
The main objective is to optimize antibody productivity and evaluate nutrient consumption under bioreactor conditions as well as minimizing medium requirements.

## Repository Contents

### **1. `CHOBiorFN_all_aa_nutrients_HP_pru1_def.py`**
- **Purpose**: Core script for simulating antibody production in CHO cells using the `iCHOv1` model.
- **Features**:
  - Defines and initializes bioreactor parameters (e.g., dilution rate, biomass concentration).
  - Adds antibody production reactions to the model using data from `added_ab_reactions.ods`.
  - Optimizes metabolic fluxes for antibody production.
  - Generates a flexible net representation (`fnet`) for simulation.
- **Key Functions**:
  - `loadCHOmodel`: Loads the metabolic model and incorporates antibody reactions.
  - `comProductivity`: Simulates antibody production and calculates maximum productivity.

### **2. `CHOBiorFN_all_aa_nutrients_HP_pru1_def_medium.py`**
- **Purpose**: Optimizes the nutrient medium composition for cost-effectiveness.
- **Features**:
  - Computes minimal nutrient uptake rates required for a given antibody production level.
  - Adds pricing constraints to simulate cost minimization for nutrient supply.
  - Incorporates bioreactor dynamics, including glucose and amino acid uptake.
- **Key Functions**:
  - `EcoMinimizeMediumim`: Optimizes the medium composition while maintaining steady-state conditions for antibody production.

### **3. `added_ab_reactions.ods`**
- **Purpose**: Provides stoichiometric coefficients for antibody production reactions.
- **Features**:
  - Contains details of reactions and their associated metabolites for antibodies like IgG.
  - Used by `loadCHOmodel` to incorporate antibody synthesis into the metabolic model.

## Workflow

1. **Prepare the Environment**:
   - Ensure all requisites are installed:
     ```bash
     pip install numpy cobra pandas matplotlib pandas-ods-reader fnyzer
     ```

2. **Load and Modify the Model**:
   - Use `CHOBiorFN_all_aa_nutrients_HP_pru1_def.py` to load the `iCHOv1` model, add antibody reactions, and simulate productivity.
   - Example:
     ```bash
     python3 CHOBiorFN_all_aa_nutrients_HP_pru1_def.py -D 0.0166 -X_ini 3.18 -X_fin 3.18 -data HP
     ```

3. **Optimize the Medium**:
   - Use `CHOBiorFN_all_aa_nutrients_HP_pru1_def_medium.py` to minimize the cost of nutrients required for antibody production.

4. **Analyze Results**:
   - The scripts output:
     - Nutrient uptake rates for optimal productivity.
     - Cost minimization results for medium composition.
     - Antibody production rates under bioreactor conditions.

## Parameters

- **Dilution Rate (`D`)**: The rate at which the medium is replaced in the bioreactor.
- **Biomass Concentration (`X_ini`, `X_fin`)**: Initial and final biomass concentrations in the reactor.
- **Dataset (`data`)**: It can take two values (HP or Late_Exp). If you want to use the dataset that contains the maximum uptake rates for glucose and amino acids in the HP experiment from [1](https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/bit.24983), use `HP`. If you want to reproduce the results from the Late Exponential Phase experiment in [2](https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/bit.24445), use `Late_Exp`.

## Solvers

Before running the simulations, it is essential to have installed a linear programming solver. In this work we used [CPLEX](https://www.ibm.com/es-es/products/ilog-cplex-optimization-studio), but others such as [GLPK](https://www.gnu.org/software/glpk/) and [Gurobi](https://www.gurobi.com/) are supported as well.

To easily install GLPK:

```bash
     sudo apt-get install python-glpk
     sudo apt-get install glpk-utils
     ```
