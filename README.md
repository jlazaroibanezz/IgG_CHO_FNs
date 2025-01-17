# Genome-Scale Metabolic Model Simulations for CHO Cells

This repository provides Python scripts to work with the genome-scale metabolic model `iCHOv1`, simulating antibody production in Chinese Hamster Ovary (CHO) cells. The [iCHOv1](http://bigg.ucsd.edu/models/iCHOv1) model is a genome-scale metabolic model of Chinese hamster ovary (CHO) cells, widely used for therapeutic protein production. It includes 1766 genes, 6663 reactions, and 4456 metabolites, covering key pathways like glycolysis, TCA cycle, and amino acid biosynthesis. The model supports optimizing antibody production, analyzing metabolic fluxes, and improving bioprocess conditions.
The main objective is to optimize antibody productivity and evaluate nutrient consumption under bioreactor conditions as well as minimizing medium requirements.

## Repository Contents

### **1. `CHOBiorFN.py`**
- **Purpose**: Core script for simulating antibody production in CHO cells using the `iCHOv1` model.
- **Features**:
  - Defines and initializes bioreactor parameters (e.g., dilution rate, biomass concentration).
  - Adds antibody production reactions to the model using data from `added_ab_reactions.ods`.
  - Optimizes metabolic fluxes for antibody production.
  - Generates a flexible net representation (`fnet`) for simulation.
- **Key Functions**:
  - `loadCHOmodel`: Loads the metabolic model and incorporates antibody reactions.
  - `comProductivity`: Simulates antibody production and calculates maximum productivity.

### **2. `CHOBiorFN_medium.py`**
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

2. **Load and Execute the Model**:
   - Use `CHOBiorFN.py` to load the `iCHOv1` model, add antibody reactions, and simulate productivity.
   - Example:
     ```bash
     python3 CHOBiorFN.py -D 0.0166 -X_ini 3.18 -X_fin 3.18 -data HP
     ```
The output of the previous command is: 

```plaintext
Growth with Xf in [ 3.18 ] gdcw L-1:
Solution: 0.00019389385950403315 mM h-1
```
3. **Optimize the Medium**:
   - Use `CHOBiorFN_medium.py` to minimize the cost of nutrients required for antibody production. This file calls the `fnet` and the `solution` generated in `CHOBiorFN.py`.  
    For example:
     ```bash
     python3 CHOBiorFN_medium.py -D 0.0166 -X_ini 3.18 -X_fin 3.18 -data HP
     ```
The output:

```plaintext
Growth with Xf in [ 3.19 ] gdcw L-1:
Solution: 0.00019389385950403315 mM h-1
Glc enters in tank:  0.0 mM
Gln enters in tank:  0.0 mM
Phe enters in tank:  1.0854803393544883 mM
Arg enters in tank:  0.0 mM
Asn enters in tank:  0.0 mM
Asp enters in tank:  0.0 mM
Cys enters in tank:  0.6385813272335126 mM
His enters in tank:  0.6103572827466555 mM
Ile enters in tank:  1.0287967078517914 mM
Leu enters in tank:  2.413129328082505 mM
Lys enters in tank:  2.3801969190054995 mM
Met enters in tank:  0.47668609675330026 mM
Pro enters in tank:  0.0 mM
Ser enters in tank:  0.0 mM
Thr enters in tank:  2.098775472431118 mM
Trp enters in tank:  0.38392230458869714 mM
Tyr enters in tank:  1.1193855421686756 mM
Val enters in tank:  2.4065777786456595 mM
Glu enters in tank:  0.0 mM
The weighted sum for the economic minimization of the medium is  0.365703463159436

```

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
