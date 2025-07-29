# -*- coding: utf-8 -*-

# load libraries
import cobra as cb
import cobra as cb
from cobra.sampling import sample
from cobra.io import load_json_model, save_json_model, write_sbml_model, read_sbml_model
import numpy as np
import random
import time
import pandas as pd
from swiglpk import *
import scanpy as sc
from scanpy import AnnData as ad
import os


# Initialize LP problem
def initilize_lp_problem(S):

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
        ia[i]=ind_row[i-1]
        ja[i]=ind_col[i-1]
        ar[i] = values[i-1]
    
    nrows=S.shape[0]
    ncol=S.shape[1]
    
    return len_vector, values, indexes, ia, ja, ar, nrows, ncol
    
    

# Solve LP problem from the structure of the metabolic model
def create_and_solve_lp_problem(lb,ub,nrows, ncol, len_vector, ia, ja, ar, 
                                obj_coefs,reactions,return_lp=None):
    
    
    lp = glp_create_prob();
    glp_set_prob_name(lp, "sample");
    glp_set_obj_dir(lp, GLP_MAX);
    glp_add_rows(lp, nrows);
    #eps=0.0000001
    eps = 1e-16
    for i in range(nrows):
        glp_set_row_name(lp, i+1, "constrain_"+str(i+1));
        glp_set_row_bnds(lp, i+1, GLP_FX, 0.0, 0.0);
    glp_add_cols(lp, ncol);
    for i in range(ncol):
        glp_set_col_name(lp, i+1, "flux_"+str(i+1));
        glp_set_col_bnds(lp, i+1, GLP_DB,lb[i]-eps,ub[i]+eps);
    glp_load_matrix(lp, len_vector, ia, ja, ar);
    
    # Solve the problem
    fluxes,Z=solve_lp_problem(lp,obj_coefs,reactions)
    
    if return_lp:
        return fluxes,Z,lp
    else:
        glp_delete_prob(lp);
        return fluxes,Z

    
# Solve LP problem from the structure of the metabolic model
def solve_lp_problem(lp,obj_coefs,reactions):
   
    # Set the coefficients of the objective function
    i=1
    for ind_coef in obj_coefs:
        glp_set_obj_coef(lp, i, ind_coef);
        i+=1

    # Initialize the parameters    
    params=glp_smcp()
    params.presolve=GLP_ON
    glp_init_smcp(params)
    params.tm_lim=120
    
    
    # Solve the problem
    glp_scale_prob(lp,GLP_SF_AUTO)
    value=glp_simplex(lp, params);
    
    #Z = glp_get_obj_val(lp);
    

    
    if value==0 or value==9:
        
        #fluxes=pd.Series(index=reactions)
        #for i in range(len(reactions)): fluxes.iloc[i]=glp_get_col_prim(lp, i+1);
        if value==9:
            raise TimeoutError("Il tempo massimo di risoluzione del problema è stato superato.")
        
        fluxes = []
        for i in range(len(reactions)): fluxes.append(glp_get_col_prim(lp, i+1))
        

        #return fluxes,Z
        return fluxes,0
    else:
        return print("error in LP problem. Problem:",str(value))


# Create LP structure
def create_lp_structure(model):
    
    reactions=[el.id for el in model.reactions]
    coefs_obj=[reaction.objective_coefficient for reaction in model.reactions]
    
    # Lower and upper bounds
    lb=[reaction.lower_bound for reaction in model.reactions]
    ub=[reaction.upper_bound for reaction in model.reactions]
    
    # Create S matrix
    S=cb.util.create_stoichiometric_matrix(model,array_type="dok")
    
    return S,lb,ub,coefs_obj,reactions

# Create random coefficients for CBS
def randomObjectiveFunction(model,executions,df_fva, start_seed=0):
    

        #reactions = model.reactions
        reactions = [reaction.id for reaction in model.reactions]
        cont=start_seed
        list_ex=reactions.copy()
        list_ex.append("type_of_problem")
        coefficients_df = pd.DataFrame(index=list_ex,columns=[str(i) for i in range(executions)])

        for i in range(0, executions):
            #df = pd.DataFrame(columns=reactions)
            #logging.info("Creating " + str(i) + "-th file - ")
            
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
                    
                    val_max=np.max([df_fva.loc[reaction,"minimum"],df_fva.loc[reaction,"maximum"]])
                    
                    if val_max!=0: #solo se la fva è diversa da zero
                        coefficients_df.loc[reaction,str(i)] = c/val_max #divido per la fva
                    else:
                        coefficients_df.loc[reaction,str(i)] = 0

                else:
                    coefficients_df.loc[reaction,str(i)] = 0

            cont=cont+1
            random.seed(cont)
                    
            if random.random()<0.5:
                coefficients_df.loc["type_of_problem",str(i)] = 0 #maximize
            else:
                coefficients_df.loc["type_of_problem",str(i)] = 1 #minimize
                            
        return coefficients_df

# CBS sampling interface to COBRA

def randomObjectiveFunctionSampling(model, nsample, coefficients_df, df_sample):

    S,lb,ub,coefs_obj,reactions = create_lp_structure(model)
    len_vector, values, indexes, ia, ja, ar, nrow, ncol = initilize_lp_problem(S)
    
    for i in range(nsample):
        #setto i coefficienti
        coefs_obj=coefficients_df.iloc[:,i].values
        if coefs_obj[-1]==1: #minimize
            coefs_obj= coefs_obj[0:-1] * -1
        else:
            coefs_obj=coefs_obj[0:-1]

    
        fluxes,Z = create_and_solve_lp_problem(lb,ub, nrow, ncol, len_vector, 
                                                        ia, ja, ar, coefs_obj,reactions,return_lp=False)
        
        df_sample.loc[i] = fluxes 

    pass

def randomObjectiveFunctionSampling_cobrapy(model, nsample, coefficients_df, df_sample, reactions):
    for i in range(nsample):
        dict_coeff={}
        if(coefficients_df.iloc[-1][i]==1):
            type_problem = -1 #minimize
        else:
            type_problem = 1
        for rxn in reactions:
            dict_coeff[model.reactions.get_by_id(rxn)] = coefficients_df.loc[rxn][i] * type_problem
        model.objective = dict_coeff
        solution =  model.optimize().fluxes
        for rxn, flux in solution.items():
            df_sample.loc[i][rxn] = flux
    pass

# CBS sampling interface to user
def model_sampling(cell, nsample, batch_size, df_coefficients, medium_cell, reactions, modelsPath, samplesPath):
    
    model2 = load_json_model(modelsPath + "/modello_cell_" + cell + ".json")
    
    if not os.path.exists(samplesPath + cell): 
        os.makedirs(samplesPath + cell)
    nBatch = int(nsample/batch_size)
    
    df_sample = pd.DataFrame(columns = reactions, index = range(batch_size))
    
    for i in range(0, nBatch):
        try:
            randomObjectiveFunctionSampling(model2, batch_size,df_coefficients.iloc[:,i*batch_size:(i+1)*batch_size], df_sample)
        except (TypeError, TimeoutError) as e:
            print(cell," " ,i, " Errore: ", e.args[0])
            randomObjectiveFunctionSampling_cobrapy(model2, batch_size,df_coefficients.iloc[:,i*batch_size:(i+1)*batch_size], df_sample, reactions)
        df_sample.to_pickle(samplesPath + cell + "/" + str(i) + ".pkl")


# CBS sampling interface to user
def model_sampling_cobrapy(cell, nsample, batch_size, df_coefficients, medium_cell, reactions, modelsPath, samplesPath):
    
    model2 = load_json_model(modelsPath + "/modello_cell_" + cell + ".json")
    
    if not os.path.exists(samplesPath + cell): 
        os.makedirs(samplesPath + cell)
    nBatch = int(nsample/batch_size)
    
    df_sample = pd.DataFrame(columns = reactions, index = range(batch_size))
    
    for i in range(0, nBatch):
        randomObjectiveFunctionSampling_cobrapy(model2, batch_size,df_coefficients.iloc[:,i*batch_size:(i+1)*batch_size], df_sample, reactions)
        df_sample.to_pickle(samplesPath + cell + "/" + str(i) + ".pkl")

# Compute fluxes averages across samples and cells
def samples_mean(cells, rxns, samplesPath):
    
    df_mean = pd.DataFrame(columns = rxns, index = cells)
    j=0
    for cell in cells:
        print(j)
        df = pd.read_pickle( samplesPath + cell + "/" + cell+ ".pkl")
        for rxn in rxns:
            df_mean.loc[cell, rxn] = df[rxn].mean()
        j+=1
    ad(df_mean.astype(float)).write_h5ad(samplesPath + 'mean')

# Compute fluxes statistics across samples and cells
def samples_stats(cells, rxns, samplesPath):
    
    columns = []
    for rxn in rxns:
        columns.append(rxn + "_min")
        columns.append(rxn + "_1q")
        columns.append(rxn + "_2q")
        columns.append(rxn + "_3q")
        columns.append(rxn + "_max")

    df_stats = pd.DataFrame(columns = columns, index = cells)
    j=0
    for cell in cells:
        print(j)
        df = pd.read_pickle( samplesPath + cell + "/" + cell+ ".pkl")
        for rxn in rxns:
            df_stats.loc[cell, rxn + "_min"] = df[rxn].min()
            df_stats.loc[cell, rxn + "_1q"] = df[rxn].quantile(0.25)
            df_stats.loc[cell, rxn + "_2q"] = df[rxn].quantile(0.50)
            df_stats.loc[cell, rxn + "_3q"] = df[rxn].quantile(0.75)
            df_stats.loc[cell, rxn + "_max"] = df[rxn].max()
        j+=1
    ad(df_stats.astype(float)).write_h5ad(samplesPath + 'stats')



