# TLOF
Transcription-based Lasso Objective Finder(TLOF) is an optimization based method to obtain a context-specific objective function for a given condition. For the complete tutorial see [here](https://github.com/opencobra/COBRA.jl/tree/master/tutorials)

This function can be called simply, by a single line of code:


```julia
using COBRA
c,obj = TLOF(metabolic_model,lambda,flux_estimation,module_flux,selected_rxns,carbon_uptake_rxn,carbon_uptake_rate,sd)
```

#### Input:
  **metabolic model**: The Metabolic model in .xml format.
  
  **lambda**: Regularization coefficient for the L1 norm term.
  
  **flux_estimation**: The flux (or estimation of flux) data which is going to be used to calcualte the context-specific objective function.
  
  It is a dataframe that has two columns, the first one contains the name of the reactions and the second one flux values.

*The next two arguments can either be given by the user or assessed by `TLOF_Preprocess` function

**module_flux**: A matrix whose dot product with the predicted flux vector returns the appropriate vector required to solve the problem.

**rxn_names**: A vector containing the name of the reactions with measured or estimated flux.

**selected_rxns**: Reactions that are meant to be included in potential cellular objective set.

**carbon_uptake_rxn**: The name of the reaction through which carbon is uptaken by a cell.  

**carbon_uptake_rate**: The exchange flux associated with the carbon source, measured experimentally.

#### OPTIONAL INPUTS

**sd**: The standard deviation value for carbon uptake rate measurement. 

  
 #### Output:

  **c**: It is the objective function found by`TLOF`.
 
  
  **obj**: The optimal value for objective function

TLOF_Preprocess can be used to find rxn_names and module_flux 


```julia
rxn_names,module_flux=TLOF_Preprocess(flux_estimation)
```
