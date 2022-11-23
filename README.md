# TLOF
Transcription-based Lasso Objective Finder(TLOF) is an optimization based method to obtain a context-specific objective function for a given condition.

## Formulation
TLOF solves the following optimization problem to find a context-specific objective function
Minimize: ‖𝒗−𝒗_𝒆𝒔𝒕 ‖<sub>𝟐</sub>+𝑹∗‖𝒄‖<sub>1</sub>
 ∑_(iϵ P)▒〖c_j v_j=uptake*g〗

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
  **metabolic model**: metabolic models conatin sotoichiometric matrix above all and also other information such as flux boundaries and Gene-Protein-Reaction rules. They can be found in different foramts including .xml. They can be downloaded from [BiGG Models](http://bigg.ucsd.edu/) or elsewhere.

  **lambda** : regularization coefficient for the L1 norm term in the objective function of the optimization problem. The larger lambda is, the more sparse the objective functions will be.
  
  **flux_estimation**: is a dataframe (or ?) that has two columns, the first one contains the name of the reactions and the second one flux values

The next two arguements can either be given by the user or assesed by `TLOF_Preprocess` function, provided in this repo.

**module_flux**:Sometimes measuring the flux of a single reactoin is not possible, thus we have measured (or estimated) flux, for example,associated with A-B or A+B where A and B are reactions in metabolic network.     fluxes are not for single reactions, but rather a module of reactions, this functions extracts the name of the reactions and also a matrix (named module flux) whose dot product whith the flux vector yields relation between single reactions and a flux module

**selected_rxns**: 

**carbon_uptake_rxn**: 

**carbon_uptake_rate**

**sd**
  
  
 #### Output:

  **c**: is of type Vector{Float64} (a vector whose elements are Float94), So this can be indexed and used like anyother vector. 
  
  **obj**: is of type Vector{Float64} (a vector whose elements are Float94), So this can be indexed and used like anyother vector.
