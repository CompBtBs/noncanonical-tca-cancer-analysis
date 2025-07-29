
#%% Libraries

from swiglpk import *
import cobra as cb
import random
import pandas as pd

#%% Create LP structure of the problem
def create_lp_structure(model):
    
    reactions=[el.id for el in model.reactions]
    coefs_obj=[reaction.objective_coefficient for reaction in model.reactions]
    
    # Lower and upper bounds
    lb=[reaction.lower_bound for reaction in model.reactions]
    ub=[reaction.upper_bound for reaction in model.reactions]
    
    # Create S matrix
    S=cb.util.create_stoichiometric_matrix(model,array_type="dok")
    
    return S,lb,ub,coefs_obj,reactions


#%% function to solve LP problem from the structure of the metabolic model
def solve_lp_problem(lp,obj_coefs,reactions):
   
    # Set the coefficients of the objective function
    i=1
    for ind_coef in obj_coefs:
        glp_set_obj_coef(lp, i, ind_coef);
        i+=1

    # Initialize the parameters    
    params=glp_smcp()
    glp_init_smcp(params)
    #params.tm_lim=10
    
    # Solve the problem
    glp_scale_prob(lp,GLP_SF_AUTO)
    value=glp_simplex(lp, params);

    if value==0 or value==9:
    
        fluxes=pd.Series(index=reactions)
        for i in range(len(reactions)): fluxes.iloc[i]=glp_get_col_prim(lp, i+1);
        
        if value==9:
            print("Time limit exceeds")
        return fluxes
    else:
        return print("error in LP problem. Problem:",str(value))
    

#%% function to solve LP problem from the structure of the metabolic model
def create_and_solve_lp_problem(S,lb,ub,obj_coefs,reactions,
                                return_lp=None):
    
    # Set the problems
    len_vector=len(S.keys())
    values=list(S.values())
    indexes=list(S.keys())
    ia = intArray(len_vector+1); 
    ja = intArray(len_vector+1);
    ar = doubleArray(len_vector+1);
    
    i=0
    ind_row=[indexes[i][0]+1 for i in range(0, len(values) )]
    ind_col=[indexes[i][1]+1 for i in range(0, len(values) )]
    for i in range(1, len(values) + 1): 
        ia[i]=int(ind_row[i-1])
        ja[i]=int(ind_col[i-1])
        ar[i] = values[i-1]
    
    #
    nrows=S.shape[0]
    ncol=S.shape[1]
    
    # 
    lp = glp_create_prob();
    glp_set_prob_name(lp, "sample");
    glp_set_obj_dir(lp, GLP_MAX);
    glp_add_rows(lp, nrows);
    eps=1e-10
    for i in range(nrows):
        glp_set_row_name(lp, i+1, "constrain_"+str(i+1));
        glp_set_row_bnds(lp, i+1, GLP_FX, 0.0, 0.0);
    glp_add_cols(lp, ncol);
    for i in range(ncol):
         glp_set_col_name(lp, i+1, "flux_"+str(i+1));
         glp_set_col_bnds(lp, i+1, GLP_DB,lb[i]-eps,ub[i]+eps);
    glp_load_matrix(lp, len_vector, ia, ja, ar);
    
    # Solve the problem
    fluxes=solve_lp_problem(lp,obj_coefs,reactions)
    
    if return_lp:
        return fluxes,lp
    else:
        glp_delete_prob(lp)
        return fluxes



def set_obj_function(model,npop, coefficients_df):
    #funzione per creare la funzione obiettivo di popolazione

    coefficients_df2=coefficients_df.iloc[:-1]
    coefficients_df2=coefficients_df2[coefficients_df2.abs()>0]
    reactions=list(coefficients_df.index)
    
    indexes=coefficients_df2.index
    reactions=[reaction for reaction in model.reactions if reaction.id.split("_cell")[0].replace("_#","") in indexes]

    type_of_problem=coefficients_df.iloc[-1]
    if type_of_problem!="maximization":
        coefficients_df2=-coefficients_df2

    dict_reactions=dict()
    for reaction in reactions:
        reaction_id=reaction.id.split("_cell")[0].replace("_#","")
        dict_reactions[reaction]=coefficients_df2.loc[reaction_id]
    model.objective=dict_reactions
    
    return model


def randomObjectiveFunction(model,executions,dfFVA,start_seed=0):
    
        dfMaxValue=dfFVA.abs().max(1)
        #reactions = model.reactions
        reactions = [reaction.id for reaction in model.reactions]
        cont=start_seed
        list_ex=reactions.copy()
        list_ex.append("type_of_problem")
        coefficients_df = pd.DataFrame(index=list_ex,columns=[str(i) for i in range(executions)])

        for i in range(0, executions):

            cont=cont+1
            random.seed(cont)
            
            # Genera un numero casuale tra 0 e 1
            threshold = random.random() #coefficiente tra 0 e 1   
    
            for reaction in reactions:
                cont=cont+1
                random.seed(cont)         
                val=random.random()   
                
                if val>threshold:
                    cont=cont+1
                    random.seed(cont)                           
                    c=2*random.random()-1 #coefficiente tra -1 e 1
                    coefficients_df.loc[reaction,str(i)] = c/dfMaxValue.loc[reaction] #divido per la fva
                else:
                    coefficients_df.loc[reaction,str(i)] = 0

            cont=cont+1
            random.seed(cont)
                    
            if random.random()<0.5:
                coefficients_df.loc["type_of_problem",str(i)] = "maximize"
            else:
                coefficients_df.loc["type_of_problem",str(i)] = "minimize"
                            
        return coefficients_df


def newObjfun(model,objective,n_pop):
    #come somma di tutte
    reactions_ids=[el.id for el in model.reactions]
    coefficients = dict()
    for reaction in model.reactions: 
        coefficients[model.reactions.get_by_id(reaction.id)] = 0
    if objective+"_#" in reactions_ids:
            coefficients[model.reactions.get_by_id(objective+"_#")] = 1
    else:
        for i in range(n_pop):
            coefficients[model.reactions.get_by_id(objective+"_cell"+str(i))] = 1                    
        
    model.objective=coefficients

    return model
