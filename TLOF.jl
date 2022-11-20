
using JuMP, Ipopt , SBML , LinearAlgebra


# flux_estimation is a dataframe (or ?) that has two columns, the first one contains the name of the reactions and the second one flux values


function TLOF(metabolic_model,lambda,flux_estimation,module_flux,selected_rxns,carbon_uptake_rxn,carbon_uptake_rate,sd=0)

    #determining irreversible reactions to set the appropriate boundary
    irreversible_indices=[]
    lb,ub=flux_bounds(metabolic_model)
    
    for i in 1:length(metabolic_model.reactions)        
        if lb[i][1]==0
            push!(irreversible_indices,i)
        end    
    end

    
    #getting the stoichiometric matrix 
    _,_,s_matrix=stoichiometry_matrix(metabolic_model)
  
   
    reaction_names_keySet=keys(metabolic_model.reactions)
    reaction_names=[]
    for item in reaction_names_keySet
        push!(reaction_names,item)

    end
    #finding the index of the carbon source uptake reaction 
    carbon_uptake_indx=findall(x->x==carbon_uptake_rxn,reaction_names )

    #Defining the optimization model and the solver 
    model=JuMP.Model(optimizer_with_attributes(Ipopt.Optimizer))
    
    #The name of the reactions are modifed to be the same as metabolic model
    modified_rxn_names=[]
    for item in rxn_names
        push!(modified_rxn_names,string("R_",item))
    end    
    #extracting the indices of reactions whose measured flux are available (this information will be used in the objective function )    
    pre_indcs_of_measured_rxns=[findall(x->x==item,reaction_names) for item in modified_rxn_names]
    indcs_of_measured_rxns=[item[1] for item in pre_indcs_of_measured_rxns ] 


    
    ommited_rxn=setdiff(1:length(metabolic_model.reactions) ,selected_rxns)   
    ommited_rxn=setdiff(ommited_rxn,carbon_uptake_indx)   

    @variable(model,-1000<=v[1:length(metabolic_model.reactions) ]<=1000)
    
    # the dual variables
    @variable(model,u[1:length(metabolic_model.species)] )
    @variable(model,g)

    # defining (the variable associated with the coefficient of each reaction in the objective function
    @variable(model,-1<=c[1:size(selected_rxns,1)]<=1)

 
    #Defining a new variable to omit absolute function from the objective function
    @variable(model,a[1:size(selected_rxns,1)])
 
 
 
    ###constraints

    @constraint(model,a.>=c)                                           
    @constraint(model,a.>=-c)    

    #Constraint regarding irreversible reactions
    @constraint(model,[j in irreversible_indices],0<=v[j]<=1000)                                     
    
    #constraint regarding the equality of the primal and the dual objective finctions, applied for reactions in potential cellular objective set
    @NLconstraint(model,sum(c[j] * v[selected_rxns[j]] for j in 1:length(selected_rxns)) ==carbon_uptake_rate*g )

    #the stoichiometric constraint, coming from the primal
    @constraint(model, s_matrix*v.==0)           

    # A constraint from the dual problem, it is applied for reactions in potential cellular objective set
    @constraint(model,[dot(u,s_matrix[:,selected_rxns[j]])-c[j] for j in 1:size(selected_rxns,1) ].>=0)        #5 ++

    # Another constraint from the dual problem, it is applied for reactions which are not in potential cellular objective set
    #@constraint(model,[dot(u,s_matrix_df[:,j+1]) for j in the_rest_indices].>=0)                                 #6 +
    @constraint(model, [j in ommited_rxn], u' * s_matrix[:, j] >= 0)                                   #6 ?

    #constraint on eschange flux for carbon source, if no standard deviation is given to the function, it will be an equality constraint
    @constraint(model,v[carbon_uptake_indx].>=[carbon_uptake_rate]-sd)                                                              
    @constraint(model,v[carbon_uptake_indx].<=[carbon_uptake_rate]-sd)                                                              


    @constraint(model,dot(u,s_matrix[:,carbon_uptake_indx[1]])+g>=0)                                               #7 +

    #the length of the c vector, which is actually the number of nonzero elements, is calculated for normalization purposes  
    n_selected=size(selected_rxns,1)

    #the length of the c vector, which is actually the number of nonzero elements, is calculated for normalization purposes  
    n_estimated=size(flux_estimation,1)

    @NLobjective(model,Min,sqrt(sum((sum(module_flux[j,i]*v[indcs_of_measured_rxns[i]] for i in 1:size(indcs_of_measured_rxns,1))- flux_estimation[j,2])^2 for j in 1:size(flux_estimation,1 )))+(n_estimated/n_selected)*lambda*(sum(a[k] for k in 1:size(c,1))))
    

    #giving random initial values to the variables
    for i in 1:size(v,1)
        if i in irreversible_indices
            set_start_value(v[i],rand(0:0.01:1000))
        else
            set_start_value(v[i],rand(-1000:0.01:1000))
        end
    end
    JuMP.optimize!(model)


    return JuMP.value.(c)
end

#bar(c)


