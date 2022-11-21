"""
When the measured (or estimated) fluxes are not for single reactions, but rather a module of reactions,
this functions extracts the name of the reactions and also a matrix (named module flux)
whose dot product whith the flux vector yields relation between single reactions and a flux module   
"""

using JuMP, Ipopt , SBML

# flux_estimation is a dataframe (or ?) that has two columns, the first one contains the name of the reactions and the second one flux values


function TLOF_Preprocess(flux_estimation)

    #the first for loop extracts the name of the reactions
    rxn_names=Vector{String31}()
    for n in 1:size(flux_estimation,1)
    
        if occursin("+" , flux_estimation[n,1]) 
            reacts=flux_estimation[n,1]
            reacts=split(reacts,"+")
            append!(rxn_names,reacts)
        elseif occursin("-" , flux_estimation[n,1])
            reacts=flux_estimation[n,1]
            reacts=split(reacts,"-")
            append!(rxn_names,reacts)

        elseif occursin("," , flux_estimation[n,1])
            reacts=flux_estimation[n,1]
            reacts=split(reacts,",")
            append!(rxn_names,reacts)
            
        else
            reacts=flux_estimation[n,1]
            push!(rxn_names,reacts)
            
        end        
    end
   
    
    #Where flux measurements are not related to a single reaction, but rather a module of reactions 
    module_flux=zeros(Float64,size(flux_estimation,1),size(rxn_names,1))

    for n in 1:size(flux_estimation,1)
    
        if occursin("+" , flux_estimation[n,1]) 
            reacts=flux_estimation[n,1]
            reacts=split(reacts,"+")
            column_indx=findall(x->x==reacts[1],rxn_names)
            column_indx2=findall(x->x==reacts[2],rxn_names)
            module_flux[n,column_indx[1]]=1
            module_flux[n,column_indx2[1]]=1
        elseif occursin("-" , flux_estimation[n,1])
            reacts=flux_estimation[n,1]
            reacts=split(reacts,"-")
            column_indx=findall(x->x==reacts[1],rxn_names)
            column_indx2=findall(x->x==reacts[2],rxn_names)
            module_flux[n,column_indx[1]]=1
            module_flux[n,column_indx2[1]]=-1

        elseif occursin("," , flux_estimation[n,1])
            reacts=flux_estimation[n,1]
            reacts=split(reacts,",")
            column_indx=[]

            #here indices of all reactions for the module are obtain but as their fluxes are considered equal, the first one is used
            for item in reacts
                push!(column_indx,findall(x->x==item,rxn_names))
            end
            
            module_flux[n,column_indx[1]]=[1]
            
        else
            reacts=flux_estimation[n,1]
            column_indx=findall(x->x==reacts,rxn_names)
            module_flux[n,column_indx]=[1]
        end
    end

    return rxn_names,module_flux 

end


