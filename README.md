# TLOF
Transcription-based Lasso Objective Finder(TLOF) is an optimization based method to obtain a context-specific objective function for a given condition.

## Formulation
TLOF solves the following optimization problem to find a context-specific objective function

Minimize:

$$ \parallel v - v<sub>\text{est}</sub> \parallel <sub>𝟐</sub> $$

‖𝒗−𝒗_𝒆𝒔𝒕 ‖<sub>𝟐</sub>+𝑹∗‖𝒄‖<sub>1</sub>

Subject to:

$$\sum_{j \in P}c_j v_j=\text{carbon uptake rate} \times  \text{g}$$

$$ a_j \geq c_j \quad \forall j \in P$$

$$ a_j \geq -c_j \quad \forall j \in P$$

$$\sum_{j \in P}S_ij v_j=0 \quad \forall i \in N$$

$$v_\text{carbon uptake rxn}=\text{carbon uptake rate}$$

$$\sum_{i=1}^N u_i S_ij \geq c_j \quad \forall j \in P$$

$$\sum_{i=1}^N u_i S_ij \geq 0 \quad \forall j \notin P , \text{carbon uptake rxn}$$

$$\sum_{i=1}^N u_i S_ij + \text{g} \geq 0 \quad \forall j \notin P , \text{carbon uptake rxn}$$


....

## Prerequisites
TLOF reads SBML models by [SBML.jl](https://github.com/LCSB-BioCore/SBML.jl), models the optimization problem by [JuMP.jl](https://github.com/jump-dev/JuMP.jl) and uses [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl) as the solver. 
[LinearAlgebra.jl](https://github.com/JuliaLang/julia/blob/master/stdlib/LinearAlgebra/src/LinearAlgebra.jl) is also required in the computations inside the function.

So these four packages are necessary to run TLOF, in addition, [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl), [CSV.jl](https://github.com/JuliaData/CSV.jl),[HTTP.jl](https://github.com/JuliaWeb/HTTP.jl) and [Test.jl](https://github.com/JuliaLang/julia/blob/master/stdlib/Test/src/Test.jl) are needed to run the test script for this function. 
They can be installed as the following example :

```
using Pkg
Pkg.add("JuMP")
```
## Usage
This function can be called as follows:
```
TLOF(metabolic_model,lambda,flux_estimation,module_flux,selected_rxns,carbon_uptake_rxn,carbon_uptake_rate,sd)
```

#### Input:
  **metabolic model**: metabolic models conatin sotoichiometric matrix above all and also other information such as flux boundaries and Gene-Protein-Reaction rules. They can be found in different formats including .xml. Metabolic models can be downloaded from [BiGG Models](http://bigg.ucsd.edu/) or elsewhere.

  **lambda** : regularization coefficient for the L1 norm term in the objective function of the optimization problem. The larger lambda is, the more sparse the objective functions will be.
  
  **flux_estimation**: is a dataframe (or ?) that has two columns, the first one contains the name of the reactions and the second one flux values

*The next two arguements can either be given by the user or assesed by `TLOF_Preprocess` function, provided in this repo.

**module_flux**:Sometimes measuring the flux of a single reactoin is not possible, thus we have measured (or estimated) flux, for example,associated with A-B or A+B where A and B are reactions in metabolic network. On the other hand, the optimization problem find flux for single reactions (in that example, A and B seperately ). But in the objective function(see foemulation section above) the difference between measured flux and the corresponding predicted value should be calculated so this `module_flux`  whose dot product whith the predicted flux vector returns the appropriate value for `V`.

**rxn_names**:This arguement is a vector containing the name of the reactions and can be different from the first column of `flux_estimation` according to the explanantions for the previous arguement.

**selected_rxns**: A user can define which reactions should be included in potential cellular objective set. This can be either all reactions of the network or any subset of the reactions, defined by their index in the stoichiometric matrix. 

**carbon_uptake_rxn**: The name of the reaction through which carbon is uptaken by a cell, for example,"R_"GLCptspp". It should matches with the raaction names of metabolic network. 

**carbon_uptake_rate** : The exchange flux associated with the carbon source, measured experimentally

**sd**: Measurements are usually performed as replicates and the average value is reported, so there is also a standard deviation value. Since problems with inequality constraints converge better, if any value is given to this arguement the capacity constraint....otherwise it will be an equlaity constraint.  
  
  
 #### Output:

  **c**: It is the objective function found by TLOF and is of type Vector{Float64} (a vector whose elements are Float94), which has the same length as the `selected_reaction`.
 
  
  **obj**: The optimal value for objective function
  
  
## TLOF_Preprocess usage
As it was explained thoroughly for `module_flux` arguement above, this function compute two input data needed to run `TLOF`: 

```rxn_names,module_flux=TLOF_Preprocess(flux_estimation)```

#### Input:
**flux_estimation**: Just the same as what was mentioned above.

 #### Output:

  **rxn_names** and **module_flux**: As explaiend earlier,what are needed for TLOF
  
