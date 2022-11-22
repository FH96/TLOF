# TLOF
Transcription-based Lasso Objective Finder


TLOF reads SBML models by [SBML.jl](https://github.com/LCSB-BioCore/SBML.jl), models the optimization problem by [JuMP.jl](https://github.com/jump-dev/JuMP.jl) and uses [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl) as the solver. 
[LinearAlgebra.jl](https://github.com/JuliaLang/julia/blob/master/stdlib/LinearAlgebra/src/LinearAlgebra.jl) is also required in the computations inside the function.

So these four packages are necessary to run TLOF, in addition HTTP.jl and Test.jl are needed to run the test script for this function. 
They can be installed as the following example :

```
using Pkg
Pkg.add("JuMP")
```
This function.. and returns the predicted flux vector

V=Ridge_FBA()
