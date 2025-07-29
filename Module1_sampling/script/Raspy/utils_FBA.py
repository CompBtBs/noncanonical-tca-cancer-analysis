# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 11:55:29 2023
Last update 21/03/2024

@author: bruno.galuzzi
"""
#%% Library
import pandas as pd
import re
import sys
from cobra.manipulation import rename_genes
from scanpy import AnnData
from cobra import  Model
import cobra as cb
import numpy as np
import copy
import random
from utils_lp import create_lp_structure,solve_lp_problem,create_and_solve_lp_problem,newObjfun,randomObjectiveFunction,set_obj_function
from ras import RAS_computation
from joblib import Parallel, delayed

"""
Function to set bounds from ras
"""
def compute_bounds_from_ras(model,rasMatrix,npop,eps):
    for i in range(npop):
        for reaction in rasMatrix.columns:
            bounds_original=model.reactions.get_by_id(reaction+"_cell"+str(i)).bounds
            valRas=rasMatrix.iloc[i].loc[reaction]

            #rimappo i ras
            bounds=(valRas*bounds_original[0]-eps,
                    valRas*bounds_original[1]+eps)
            
            model.reactions.get_by_id(reaction+"_cell"+str(i)).bounds=bounds
        
    return model
""" 
Function to normalize ras and compute ras matrix
"""
def compute_ras_matrix(ras_adata,type_ras_normalization):
    #ras normalization
    Matrix=ras_adata.to_df()   
    if type_ras_normalization=="sum":
        Matrix=Matrix.div(Matrix.sum())
    else:
        Matrix=Matrix.div(Matrix.max())
    Matrix=Matrix.fillna(0)
    
    return Matrix
"""
Function to modify exchange for a population of models
"""    
    
def modify_exchange(model,exchanges,npop):
    model_modify=model.copy()
    
    #modify exchanges
    for reaction in exchanges: 
        bounds=model_modify.reactions.get_by_id(reaction).bounds
        bounds=(bounds[0]*npop,bounds[1]*npop)
        model_modify.reactions.get_by_id(reaction).bounds=bounds
    
    return model_modify
"""
Create dataframe in which reactions are rows and cells are columns

"""
        
def transformDf(df,npop):
    
    index=df.index
    ex=[el.split("_cell")[0] for el in index]
    reactions=list(set(index))

    dfCells=pd.DataFrame(index=reactions,columns=["cell"+str(i) for i in range(npop)])


    for reaction in reactions:
        if "_#" in reaction:
            dfCells.loc[reaction,:]=df.loc[reaction]
        else:
            valori=[df.loc[reaction+"_cell"+str(i)] for i in range(npop)]
            dfCells.loc[reaction,:]=valori

    return dfCells

""" 
function to create a population of models
"""

def popModel(model_orig,
             n_pop,                       #how many models to create
             objective=None,              #objective functions
             c_round=10,
             fraction_of_optimum=0,
             compute_fva=True,
             npop_fva=2,
             processes=1):

    #copy the model
    model=model_orig.copy()

    #create a model
    model_merged=createPopModel(model,n_pop)

    # create objective function
    model_merged=newObjfun(model_merged,objective,n_pop)
    

    if compute_fva:
        
        #create a model

        model_merged2=createPopModel(model_orig.copy(),npop_fva)
        
        # create objective function
        model_merged2=newObjfun(model_merged2,objective,npop_fva)
        

        
        reactions=[el.id for el in model.reactions]
        reactions_ex=[el.id for el in model.exchanges]
        reaction_noex=[el for el in reactions if el not in reactions_ex]
        reaction_list=[el+"_cell0" for el in reaction_noex]
        reaction_list.extend([el+"_#" for el in reactions_ex])
        
        dfFVA=cb.flux_analysis.flux_variability_analysis(model_merged2,
                                                         fraction_of_optimum=fraction_of_optimum,
                                                         processes=processes,
                                                         reaction_list=reaction_list                           
                                                                                    ).round(c_round)
    
        for reaction in reaction_noex:

            for i in range(n_pop):
                  model_merged.reactions.get_by_id(reaction+"_cell"+str(i)).bounds=(dfFVA.loc[reaction+"_cell0","minimum"],
                                                                                   dfFVA.loc[reaction+"_cell0","maximum"])
        for reaction in reactions_ex:
                  model_merged.reactions.get_by_id(reaction+"_#").bounds=(dfFVA.loc[reaction+"_#","minimum"],
                                                                          dfFVA.loc[reaction+"_#","maximum"])                                                                 
        return model_merged, dfFVA
    else:
        return model_merged


def createPopModel(model,n_pop):
    #
    # Create a population of networks
    #
    #
    exchanges=[el.id for el in model.exchanges]
    reactions_noexchange=[el for el in model.reactions if el.id not in exchanges]
    reactions_exchange=[el for el in model.reactions if el.id  in exchanges]
    
    #considero i metaboliti delle exchanges
    metabolite_exchange=[]
    for exchange in reactions_exchange:
        for metabolite in exchange.metabolites:
            metabolite_exchange.append(metabolite.id)
    metabolite_exchange=list(set(metabolite_exchange))
    
    for reaction in reactions_exchange:
        reaction.id=reaction.id+"_#"              #gli do un nome a parte
    
    #remove exchange reactions from the model
    model.remove_reactions(reactions_exchange)  
        
    model_merged=Model("model_merged")
    reactions_to_add=[]

    for i in range(n_pop):
        print(i)
        reactions2=copy.deepcopy(reactions_noexchange)
        cell_name="_cell"+str(i)
        for reaction in reactions2:
            reaction.id+=cell_name
            for el in reaction.metabolites: #change also the metabolite names
                if "_cell" not in el.id and el.id not in metabolite_exchange:
                    el.id+=cell_name
                
        reactions_to_add.extend(reactions2)
    
    #add reactions of the super-network
    model_merged.add_reactions(reactions_to_add)
    
    #add the exchange
    model_merged.add_reactions(reactions_exchange)

    return model_merged


"""
Function to compute scFBA with cooperation 
"""
def FBA_coop(
          model_orig,                         #original model
          npop,                               #how many elements of the population
          ras_adata,                          #ras to include in the popultion for constraints
          compute_fva=True,
          type_ras_normalization="sum",
          eps=0,
          return_adata=True
          ):
    
    #single FBA con modello unico
    model=model_orig.copy()

    #compute ras matrix
    rasMatrix=compute_ras_matrix(ras_adata,type_ras_normalization)

    #constraints for each cell and reaction
    model=compute_bounds_from_ras(model,rasMatrix,npop,eps)
    
    #optimize model       
    S,lb,ub,coefs_obj,reactions=create_lp_structure(model)

    print("compute optimization:")
    fluxes=create_and_solve_lp_problem(S,lb,ub,coefs_obj,reactions,return_lp=False)
    
    #obj_res=model.optimize()  
    dfTot=transformDf(fluxes.round(10),npop)
    
    if return_adata:    
        adata=AnnData(dfTot.T)
        return adata
    else:
        return dfTot


"""
Function to compute scFBA without cooperation
"""
def FBA_nocoop(model_orig,
          objective,
          npop,
          ras_adata,
          type_ras_normalization="sum",
          compute_fva=True,
          eps=0,
          round_c=10,
          processes=1,
          fraction_of_optimum=0,
          return_adata=True,
          print_cell_computation=False
          ):
    
    # compute ras matrix
    rasMatrix=compute_ras_matrix(ras_adata,type_ras_normalization)   
    
    # Set up the bounds using FVA
    model_orig2=model_orig.copy()
    reactions=[el.id for el in model_orig2.reactions]
    
    if compute_fva:
        dfFVA=cb.flux_analysis.flux_variability_analysis(model_orig2,fraction_of_optimum=fraction_of_optimum,processes=processes).round(round_c)
        for reaction in reactions:
            model_orig2.reactions.get_by_id(reaction).bounds=(dfFVA.loc[reaction,"minimum"],dfFVA.loc[reaction,"maximum"])
    

    # Compute FBA for any sample
    dfTot=pd.DataFrame(index=[reaction.id for reaction in model_orig.reactions],
                       columns=["cell"+str(i) for i in range(npop)])

    for i in range(npop):
        if print_cell_computation and i % 10 == 0:
            print(i)
        model=model_orig2.copy()
        ras_sample=rasMatrix.iloc[i]
        #set ras values as bounds
        for reaction in rasMatrix.columns:
            bounds_original=model.reactions.get_by_id(reaction).bounds
            valRas=ras_sample.loc[reaction]

            bounds=(valRas*bounds_original[0]-eps,
                    valRas*bounds_original[1]+eps)
            
            model.reactions.get_by_id(reaction).bounds=bounds
            
            
        #optimize model       
        S,lb,ub,coefs_obj,reactions=create_lp_structure(model)
        fluxes=create_and_solve_lp_problem(S,lb,ub,coefs_obj,reactions,return_lp=False)
        dfTot.loc[:,"cell"+str(i)]=fluxes.round(round_c)
            
    #choose data type
    if return_adata:
        adata=AnnData(dfTot.T)
        return adata
    else:
        return dfTot


"""
Function to compute scFBA without cooperation RETURN models and not fluxes
"""
def FBA_nocoop_save_models(model_orig,
          objective,
          npop,
          ras_adata,
          type_ras_normalization="sum",
          compute_fva=True,
          eps=0,
          round_c=10,
          processes=1,
          fraction_of_optimum=0,
          return_adata=True,
          print_cell_computation=False,
          path_models = ""
          ):
    
    # compute ras matrix
    rasMatrix=compute_ras_matrix(ras_adata,type_ras_normalization)   
    
    # Set up the bounds using FVA
    model_orig2=model_orig.copy()
    reactions=[el.id for el in model_orig2.reactions]
    
    if compute_fva:
        dfFVA=cb.flux_analysis.flux_variability_analysis(model_orig2,fraction_of_optimum=fraction_of_optimum,processes=processes).round(round_c)
        dfFVA.to_csv(f"{path_models}/dfFVA.csv")
        for reaction in reactions:
            model_orig2.reactions.get_by_id(reaction).bounds=(dfFVA.loc[reaction,"minimum"],dfFVA.loc[reaction,"maximum"])
    
    
    # Compute FBA for any sample
    dfTot=pd.DataFrame(index=[reaction.id for reaction in model_orig.reactions],
                       columns=["cell"+str(i) for i in range(npop)])

    for i, cell_id in enumerate(rasMatrix.index):
        if print_cell_computation and i % 10 == 0:
            print(i)
        model=model_orig2.copy()
        ras_sample=rasMatrix.loc[cell_id]
        #set ras values as bounds
        for reaction in rasMatrix.columns:
            bounds_original=model.reactions.get_by_id(reaction).bounds
            valRas=ras_sample.loc[reaction]

            bounds=(valRas*bounds_original[0]-eps,
                    valRas*bounds_original[1]+eps)
            
            model.reactions.get_by_id(reaction).bounds=bounds
            
            
        #save model       
        cb.io.save_json_model(model, f"{path_models}/modello_cell_{cell_id}.json")  # JSON

    print("Models saved")

#"""
#Function to compute scFBA without cooperation RETURN models and not fluxes multi medium
#"""
#def FBA_nocoop_save_models_multi_medium(model_orig,
#          dfMedium,
#          objective,
#          npop,
#          ras_adata,
#          batch_list,
#          type_ras_normalization="sum",
#          compute_fva=True,
#          eps=0,
#          round_c=10,
#          processes=1,
#          fraction_of_optimum=0,
#          return_adata=True,
#          print_cell_computation=False,
#          path_models = ""
#          ):
#    
#    print("multi medium")
#    
#    # compute ras matrix
#    rasMatrix=compute_ras_matrix(ras_adata,type_ras_normalization)#

#    if compute_fva:
#        dfFVA=cb.flux_analysis.flux_variability_analysis(model_orig,fraction_of_optimum=fraction_of_optimum,processes=processes).round(round_c)
#        dfFVA.to_csv(f"{path_models}/dfFVA.csv")

#    for i, cell_id in enumerate(rasMatrix.index):
#
#       # Set up the bounds using FVA
#        model_orig2=model_orig.copy()

#        if batch_list[cell_id] == '2':
#            model_orig2.medium = dfMedium.to_dict()['mIVC1:adv_DMEM:F12 = 1:1']
#        else:
#            model_orig2.medium = dfMedium.to_dict()['E8_ess']

#        reactions=[el.id for el in model_orig2.reactions]
    
#        if compute_fva:
#            dfFVA=cb.flux_analysis.flux_variability_analysis(model_orig2,fraction_of_optimum=fraction_of_optimum,processes=processes).round(round_c)
#            for reaction in reactions:
#                model_orig2.reactions.get_by_id(reaction).bounds=(dfFVA.loc[reaction,"minimum"],dfFVA.loc[reaction,"maximum"])

            
        # Compute FBA for any sample
#        dfTot=pd.DataFrame(index=[reaction.id for reaction in model_orig.reactions],
#                        columns=["cell"+str(i) for i in range(npop)])

#        if print_cell_computation and i % 10 == 0:
#            print(i)
#        model=model_orig2.copy()
#        ras_sample=rasMatrix.loc[cell_id]
#        #set ras values as bounds
#        for reaction in rasMatrix.columns:
#            bounds_original=model.reactions.get_by_id(reaction).bounds
#            valRas=ras_sample.loc[reaction]
#
#            bounds=(valRas*bounds_original[0]-eps,
#                    valRas*bounds_original[1]+eps)
#            model.reactions.get_by_id(reaction).bounds=bounds
            
            
        #save model       
#        cb.io.save_json_model(model, f"{path_models}/modello_cell_{cell_id}.json")  # JSON

#    print("Models saved")




def FBA_nocoop_save_models_multi_medium(model_orig,
          dfMedium,
          objective,
          npop,
          ras_adata,
          batch_list,
          type_ras_normalization="sum",
          compute_fva=True,
          eps=0,
          round_c=10,
          processes=100,
          fraction_of_optimum=0,
          return_adata=True,
          print_cell_computation=False,
          path_models=""
          ):
    
    def process_cell(cell_id, ras_sample, model_orig, dfFVA, eps, path_models, batch):
        model = model_orig.copy()
        
        # Set medium
        if batch == '2':
            model.medium = dfMedium.to_dict()['mIVC1:adv_DMEM:F12 = 1:1']
        else:
            model.medium = dfMedium.to_dict()['E8_ess']
        
        # Set FVA bounds
        for reaction_id in dfFVA.index:
            model.reactions.get_by_id(reaction_id).bounds = (dfFVA.loc[reaction_id, "minimum"], dfFVA.loc[reaction_id, "maximum"])

        # Set ras values as bounds
        for reaction_id in ras_sample.index:
            bounds_original = model.reactions.get_by_id(reaction_id).bounds
            valRas = ras_sample[reaction_id]
            bounds = (valRas * bounds_original[0] - eps, valRas * bounds_original[1] + eps)
            model.reactions.get_by_id(reaction_id).bounds = bounds

        # Save model
        cb.io.save_json_model(model, f"{path_models}/modello_cell_{cell_id}.json")

    print("multi medium")
    
    # Compute ras matrix
    rasMatrix = compute_ras_matrix(ras_adata, type_ras_normalization)

    # Pre-calculate FVA for each medium
    if compute_fva:
        fva_results = {}
        
        # FVA for the original model
        dfFVA_orig = cb.flux_analysis.flux_variability_analysis(model_orig, fraction_of_optimum=fraction_of_optimum, processes=1).round(round_c)
        dfFVA_orig.to_csv(f"{path_models}/dfFVA_orig.csv")
        fva_results['orig'] = dfFVA_orig
        
        # FVA for the model with medium mIVC1:adv_DMEM:F12 = 1:1
        model_mIVC1 = model_orig.copy()
        model_mIVC1.medium = dfMedium.to_dict()['mIVC1:adv_DMEM:F12 = 1:1']
        dfFVA_mIVC1 = cb.flux_analysis.flux_variability_analysis(model_mIVC1, fraction_of_optimum=fraction_of_optimum, processes=1).round(round_c)
        dfFVA_mIVC1.to_csv(f"{path_models}/dfFVA_mIVC1.csv")
        fva_results['mIVC1'] = dfFVA_mIVC1

        # FVA for the model with medium E8_ess
        model_E8_ess = model_orig.copy()
        model_E8_ess.medium = dfMedium.to_dict()['E8_ess']
        dfFVA_E8_ess = cb.flux_analysis.flux_variability_analysis(model_E8_ess, fraction_of_optimum=fraction_of_optimum, processes=1).round(round_c)
        dfFVA_E8_ess.to_csv(f"{path_models}/dfFVA_E8_ess.csv")
        fva_results['E8_ess'] = dfFVA_E8_ess

    # Parallel processing of cells
    results = Parallel(n_jobs=processes)(
        delayed(process_cell)(
            cell_id,
            rasMatrix.loc[cell_id],
            model_orig,
            fva_results['mIVC1'] if batch_list[cell_id] == '2' else fva_results['E8_ess'],
            eps,
            path_models,
            batch_list[cell_id]
        ) for cell_id in rasMatrix.index
    )

    print("Models saved")


"""
Function to compute sampling with cooperation
"""
def coop_sampling(model_orig,
          npop,
	      coefficients_df,
          ras_adata,
          dfFVA,
          sample_size=100,
          type_ras_normalization="sum",
          eps=0,
          round_c=10,
          return_adata=True
          ):
    
    #single FBA con modello unico
    model=model_orig.copy()

    #ras normalization
    rasMatrix=compute_ras_matrix(ras_adata,type_ras_normalization)
    
    #constraints for each cell and reaction
    compute_bounds_from_ras(model,rasMatrix,npop,eps)
            
    sample_size=coefficients_df.shape[1]      
    
    S,lb,ub,coefs_obj,reactions=create_lp_structure(model)
    
    #create df for sampling
    df=pd.Series(index=reactions,data=0)
    
    for i in range(sample_size):
        print("compute iteration:",i)
        coefs_obj=coefficients_df.iloc[:,i].values
        if coefs_obj[-1]=="minimize":
            coefs_obj=-coefs_obj[0:-1]
        else:
            coefs_obj=coefs_obj[0:-1]
            
        #if i==0:
        fluxes=create_and_solve_lp_problem(S,lb,ub,coefs_obj,reactions,return_lp=False)
        #else:
        #    fluxes,_=solve_lp_problem(lp_problem,coefs_obj,reactions)

        df=df+fluxes.round(round_c)
        
    df=df/sample_size
    
    print("final assemble")
    dfTot=transformDf(df,npop)
    
    if return_adata:    
        adata=AnnData(dfTot.T)
        return adata
    else:
        return dfTot

"""
Main function
"""
def popFBA(model_orig,                      #metabolic network
            adata,                           #Ann data object
            objective=None,
            type_of_problem="maximization",
            nsample=100,
            cooperation=False,
            compute_fva=False,
            eps=0,
            npop_fva=2,
            type_of_data="transcriptomics",  #type of transcroptimics data (ras or rna-seq)
            type_ras_normalization="max",
            and_expression=np.nanmin,
            or_expression=np.nansum,
            fraction_of_optimum=0,
            processes=1,
            print_cell_computation=False,
            round_c=10,
            ):
    
    #Check for possible input error
    if type_of_problem not in ["maximization","sampling"]:
        return print("Error: the problem must be maximization or sampling")
    
    if type_of_data not in ["ras","transcriptomics"]:
        return print("Error: the data must be transcriptomics or ras")    

    if objective ==None and type_of_problem=="maximization":
        return print("Specify the objective function")
        #sampling
    if nsample==None and type_of_problem=="sampling":
        return print("You need to specify the sample size")
    
    if type(adata)!=type(AnnData()) and type(adata)!=type(pd.DataFrame()):
        return print("You need to AnnData or Pandas DataFrame!")
    else:
        if type(adata)==type(pd.DataFrame()):
            adata=AnnData(adata)
    
    npop=adata.shape[0]
    model=model_orig.copy()   
    exchanges=[el.id for el in model.exchanges]

    if type_of_data=="transcriptomics":
        # Compute RAS for any sample
        print("RAS computation")
        ras_object=RAS_computation(adata,model)
        ras_adata=ras_object.compute(
                                      and_expression=and_expression,
                                      or_expression=or_expression
                                     )
    else:
        ras_adata=adata.copy()
    
    #Compute fluxes 
    if type_of_problem=="maximization":
        #maximization
        if cooperation:
            #cooperation between cells
            
            #modify exchange reaction
            model_x_pop=modify_exchange(model,exchanges,npop)
            
            #create super-network of npop
            model_merged,dfFVA=popModel(model_x_pop,
                                  npop,
                                  objective,
                                  compute_fva=compute_fva,
                                  npop_fva=npop_fva
                                  )
    
            #compute fluxes
            adata_fluxes_pop=FBA_coop(model_merged,
                                npop,
                                ras_adata,
                                eps=eps,
                                type_ras_normalization=type_ras_normalization,
                                return_adata=True)
        else:
            #no cooperation between cells
            adata_fluxes_pop=FBA_nocoop(model,
                              objective,
                              npop,
                              ras_adata,
                              eps=eps,
                              compute_fva=compute_fva,
                              type_ras_normalization=type_ras_normalization,
                              return_adata=True,
                              fraction_of_optimum=fraction_of_optimum,
                              print_cell_computation=print_cell_computation,
                              processes=processes,
                              round_c=round_c
                              )
      
    else:
            
         if cooperation:
             #modify exchange
             model_x_pop=modify_exchange(model,exchanges,npop)
             
             #create super-network
             model_merged,dfFVA=popModel(model_x_pop,
                                   npop,
                                   objective,
                                   compute_fva=compute_fva,
                                   npop_fva=npop_fva
                                   )
             #create random objective coefficients
             print("Create random objective functions")
             coefficients_df=randomObjectiveFunction(model_x_pop,nsample,dfFVA,0)

             #compute fluxes
             print("Compute sampling of random objective function")
             adata_fluxes_pop=coop_sampling(model_merged,
                                 npop,
        				         coefficients_df,
                                 ras_adata,
                                 dfFVA,
                                 eps=eps,
                                 type_ras_normalization=type_ras_normalization,
                                 return_adata=True)
         else:
             #no cooperation between cells
             return print("to do")

    
    return adata_fluxes_pop


"""
Main function return models
"""
def popFBA_models(model_orig,                      #metabolic network
            adata,                           #Ann data object
            dfMedium,
            objective=None,
            type_of_problem="maximization",
            nsample=100,
            cooperation=False,
            compute_fva=False,
            eps=0,
            npop_fva=2,
            type_of_data="transcriptomics",  #type of transcroptimics data (ras or rna-seq)
            type_ras_normalization="max",
            and_expression=np.nanmin,
            or_expression=np.nansum,
            fraction_of_optimum=0,
            processes=1,
            print_cell_computation=False,
            round_c=10,
            path_models="",
            multi_med = False
            ):
    
    #Check for possible input error
    if type_of_problem not in ["maximization","sampling"]:
        return print("Error: the problem must be maximization or sampling")
    
    if type_of_data not in ["ras","transcriptomics"]:
        return print("Error: the data must be transcriptomics or ras")    

    if objective ==None and type_of_problem=="maximization":
        return print("Specify the objective function")
        #sampling
    if nsample==None and type_of_problem=="sampling":
        return print("You need to specify the sample size")
    
    if type(adata)!=type(AnnData()) and type(adata)!=type(pd.DataFrame()):
        return print("You need to AnnData or Pandas DataFrame!")
    else:
        if type(adata)==type(pd.DataFrame()):
            adata=AnnData(adata)
    
    npop=adata.shape[0]
    model=model_orig.copy()   
    exchanges=[el.id for el in model.exchanges]

    if type_of_data=="transcriptomics":
        # Compute RAS for any sample
        print("RAS computation")
        ras_object=RAS_computation(adata,model)
        ras_adata=ras_object.compute(
                                      and_expression=and_expression,
                                      or_expression=or_expression
                                     )
    else:
        ras_adata=adata.copy()
    ras_adata.write_h5ad(f"{path_models}/ras_adata.h5ad",compression='gzip')
    #Compute fluxes 
    if type_of_problem=="maximization":
        #maximization
        if cooperation:
            #cooperation between cells
            
            #modify exchange reaction
            model_x_pop=modify_exchange(model,exchanges,npop)
            
            #create super-network of npop
            model_merged,dfFVA=popModel(model_x_pop,
                                  npop,
                                  objective,
                                  compute_fva=compute_fva,
                                  npop_fva=npop_fva
                                  )
    
            #compute fluxes
            adata_fluxes_pop=FBA_coop(model_merged,
                                npop,
                                ras_adata,
                                eps=eps,
                                type_ras_normalization=type_ras_normalization,
                                return_adata=True)
        else:
            if not multi_med:
                FBA_nocoop_save_models(model,
                                objective,
                                npop,
                                ras_adata,
                                eps=eps,
                                compute_fva=compute_fva,
                                type_ras_normalization=type_ras_normalization,
                                return_adata=True,
                                fraction_of_optimum=fraction_of_optimum,
                                print_cell_computation=print_cell_computation,
                                processes=processes,
                                round_c=round_c,
                                path_models=path_models
                                )
            else:
                FBA_nocoop_save_models_multi_medium(model,
                                dfMedium,
                                objective,
                                npop,
                                ras_adata,
                                batch_list = adata.obs['batch'],
                                eps=eps,
                                compute_fva=compute_fva,
                                type_ras_normalization=type_ras_normalization,
                                return_adata=True,
                                fraction_of_optimum=fraction_of_optimum,
                                print_cell_computation=print_cell_computation,
                                processes=processes,
                                round_c=round_c,
                                path_models=path_models
                                )


            return "ok"       
    else:
            
         if cooperation:
             #modify exchange
             model_x_pop=modify_exchange(model,exchanges,npop)
             
             #create super-network
             model_merged,dfFVA=popModel(model_x_pop,
                                   npop,
                                   objective,
                                   compute_fva=compute_fva,
                                   npop_fva=npop_fva
                                   )
             #create random objective coefficients
             print("Create random objective functions")
             coefficients_df=randomObjectiveFunction(model_x_pop,nsample,dfFVA,0)

             #compute fluxes
             print("Compute sampling of random objective function")
             adata_fluxes_pop=coop_sampling(model_merged,
                                 npop,
        				         coefficients_df,
                                 ras_adata,
                                 dfFVA,
                                 eps=eps,
                                 type_ras_normalization=type_ras_normalization,
                                 return_adata=True)
         else:
             #no cooperation between cells
             return print("to do")

    
    return adata_fluxes_pop


def convert_genes(model,annotation):
    model2=model.copy()
    try:
        dict_genes={gene.id:gene.notes[annotation]  for gene in model2.genes}
    except:
        print("No annotation in gene dict!")
        return -1
    rename_genes(model2,dict_genes)

    return model2