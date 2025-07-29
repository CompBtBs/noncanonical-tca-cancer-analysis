# -*- coding: utf-8 -*-
"""
Created on Sat Jun 17 15:51:50 2023

@author: bruno.galuzzi

"""
#%% load libraries
import cobra as cb
import numpy as np
import scanpy as sc
import pandas as pd
import sys
from utils_FBA import popFBA,convert_genes
from scanpy import AnnData as ad
import magic

#%% load cobra model
# G:\.shortcut-targets-by-id\1wD75Mjo9VYOD03yLCXuo7NKUcXMrL9U5\compBTBS\Modelli\ENGRO2\ENGRO2_21032024.xml
model_orig = cb.io.read_sbml_model("Modelli/ENGRO2/ENGRO2_17042024.xml")

#%% translate gene annotation
model=convert_genes(model_orig,
                    "ENSG")

#%% load single cell RNA-seq dataset
# C:\Users\f.lapi\Documents\GitHub\scFBApy\datasets "Dati/popFBA/BC04dataset
adata=sc.read_h5ad("Dati/AnalisiDatiPadova/GSE136447_555_cln.h5ad")

#%% denoising with MAGIC agorithm
sc.external.pp.magic(adata,knn=3)
#adata=adata[1:,:]

#%% Compute optimal fluxes
adata_fluxes_pop=popFBA(model,adata,
                        cooperation=False,
                        type_of_problem="maximization",
                        objective="Biomass",
                        print_cell_computation=True,
                        path_models = ""
                        )


 #%%salvo il file
results_file_adata = "prova.h5ad" 
results_file_csv = "prova.csv"

#adata_fluxes_pop = ad(adata_fluxes_pop.to_df().round(6).T.drop_duplicates().T)
#adata_fluxes_pop.write(results_file_adata)

df_fluxes_pop = adata_fluxes_pop.to_df()
df_fluxes_pop = df_fluxes_pop.round(6)
df_fluxes_pop.to_csv(results_file_csv)

## TODO ##
# controllo qualità dati RNA
# ORA MEDIUM TUTTO APERTO -> MEDIUM DI CRESCITA DELLE CELLULE SE SI SA LA COMPOSIZIONE, ALTRIMENTI CAPIRE IL MEDIUM PIU' SENSATO
# DENOISER PER I DATI
## SPARSITA' DEL DATO SINGLE CELL
## VERIFICARE CHE NON CI SIANO CELLULE CHE CRESCONO CON I RAS SENZA DENOISER
# riprodurre i loro cluster per vedere se famo uguali