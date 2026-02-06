# Mechanistically Informed Machine Learning Links Non-Canonical TCA Cycle Activity to Warburg Metabolism and Hallmarks of Malignancy

This repository contains the code and relevant resources used in the analysis presented in the paper:

Lin L, Lapi F, Galuzzi BG, Vanoni M, Alberghina L, Damiani C (2025) Mechanistically informed machine learning links non-canonical TCA cycle activity to Warburg metabolism and hallmarks of malignancy. PLoS Comput Biol 21(12): e1013384. https://doi.org/10.1371/journal.pcbi.1013384

## Paper abstract

Cancer cells undergo extensive metabolic rewiring to support growth, survival, and phenotypic plasticity. A non-canonical variant of the tricarboxylic acid (TCA) cycle, characterized by mitochondrial-to-cytosolic citrate export, has emerged as critical for embryonic stem cell differentiation. However, its role in cancer remains poorly understood.

Here, we present a two-step computational framework to systematically analyze the activity of this non-canonical TCA cycle across over 500 cancer cell lines and investigate its role in shaping hallmarks of malignancy. First, we applied constraint-based modeling to infer cycle activity, defining two complementary metrics: *Cycle Propensity*, measuring the likelihood of its engagement in each cell line, and *Cycle Flux Intensity*, quantifying average flux through the reaction identified as rate-limiting. We identified distinct tumor-specific patterns of pathway utilization. Notably, cells with high *Cycle Propensity* preferentially rerouted cytosolic citrate via aconitase 1 (ACO1) and isocitrate dehydrogenase 1 (IDH1), promoting α-ketoglutarate (αKG) and NADPH production. Elevated engagement of this cycle strongly correlated with Warburg-like metabolic shifts, including decreased oxygen consumption and increased lactate secretion.

In the second step, to uncover non-metabolic transcriptional signatures associated with non-canonical TCA cycle activity, we performed machine learning–based feature selection using ElasticNet and XGBoost, identifying robust gene signatures predictive of cycle activity. Over-representation analysis revealed enrichment in genes involved in metastatic behavior, angiogenesis, stemness, and key oncogenic signaling. SHapley Additive exPlanations (SHAP) further prioritized genes with the strongest predictive contributions, highlighting candidates for experimental validation. Correlation analysis of DepMap gene-dependency profiles revealed distinct vulnerability patterns associated with non-canonical TCA cycle activity, outlining a characteristic landscape of genetic dependencies.

Together, our integrative framework uniting constraint-based metabolic modeling and machine learning systematically reveals how non-canonical TCA cycle dynamics underpin metabolic plasticity and promote malignant traits.

## Author Summary

A non-canonical variant of the TCA cycle that bypasses several steps of the mitochondrial TCA cycle has recently been identified. This pathway involves the export of mitochondrial citrate to the cytosol, its conversion to oxaloacetate and subsequently to malate, followed by the re-import of malate into mitochondria to complete the cycle. This non-canonical mitochondrial–cytosolic cycle has been linked to changes in stem cell identity. However, its functional role in tumor metabolism remains poorly understood.

To explore this connection, we developed a computational framework that combines mechanistic metabolic modeling with machine learning to analyze the activity of this pathway across more than 500 cancer cell lines. Using this approach, we quantified the activity of the non-canonical TCA cycle and identified gene expression programs that predict its engagement. We found that cells utilizing this pathway exhibit a Warburg-like metabolic profile. The associated gene expression programs are enriched for processes related to core hallmarks of malignancy, including metastatic behavior, angiogenesis, stemness, and key oncogenic signaling. Overall, our results demonstrate that this recently discovered pathway links metabolic rewiring with transcriptional programs that drive tumor aggressiveness and progression, suggesting a mechanistic connection between the activation of the non-canonical TCA cycle and transcriptional states associated with malignancy.

## Pipeline overview

<p align="center">
  <img src="Graphical_Abstract.png" alt="Graphical_Abstract" width="700"/>
</p>

We conceived a two-step computational strategy integrating a metabolic modeling module and a machine learning module to dissect the activity of the non-canonical TCA cycle (Cit-Mal) across cancer cell lines.

In the **metabolic modeling module**, we employed state-of-the-art constraint-based modeling techniques to infer Cit-Mal Cycle activity directly from gene expression data. Cellular metabolism is represented as a network of biochemical reactions subject to stoichiometric mass balance and physicochemical constraints, such as minimal and maximal reaction fluxes. These flux boundaries can be personalized using omics data—such as transcriptomic profiles—to reflect the molecular context of individual samples.

In this study, we tailored the boundaries of the generic core metabolic network ENGRO2.2 based on the transcriptomes of 513 cancer cell lines cultured under standardized conditions. Specifically, we computed Reaction Activity Scores (RAS) from gene expression using Gene–Protein–Reaction (GPR) rules, and used these scores to modulate the upper bounds of each reaction’s feasible flux. Each reconstructed model defines the space of feasible steady-state flux distributions for a specific cell line.

Rather than selecting a single flux vector by optimizing a predefined biological objective (e.g., biomass yield), we adopted flux sampling to systematically explore the entire solution space consistent with network structure and molecular constraints. We used a corner-based sampling algorithm to generate thousands of feasible flux distributions per model, thereby approximating the distribution of steady-state flux configurations for each cell line.

From these sampled flux distributions, we devised two complementary metrics to quantify Cit-Mal Cycle activity:

- `Cycle Propensity`: the fraction of sampled steady-state flux distributions in which all three hallmark reactions (SLC25A1, ACLY, MDH1) are simultaneously active. This metric is dimensionless (0–1) and reflects the likelihood that the cycle operates under the given constraints.
- `Cycle Flux Intensity`: the average flux through the bottleneck reaction (typically ACLY) in those sampled states where the cycle is active. This metric is expressed in arbitrary flux units, reflecting the relative nature of transcriptomics-derived constraints.

In the **machine learning module**, these *in silico*–derived metrics serve as training targets for supervised regression models designed to predict `Cycle Propensity` and `Cycle Flux Intensity` from non-metabolic gene expression data. After model training, we performed feature selection to identify robust transcriptional predictors of Cit-Mal Cycle activity, thereby revealing broader gene-expression programs potentially regulating or co-occurring with its engagement. Finally, we applied SHAP analysis to quantify the contribution of individual genes to model predictions, providing a biologically interpretable link between transcriptional features and Cit-Mal Cycle activity.



## Data availability

* **Gene expression data** were obtained from the [Cancer Cell Line Encyclopedia (CCLE)](https://sites.broadinstitute.org/ccle/datasets).

* **Validation fluxes**  
  For the validation of our cycle activity metrics, we used fluxes generated in Maspero *et al.* (2024), derived from spatial transcriptomics data from clear cell renal cell carcinoma tissue and adjacent normal parenchyma (sample ID: `PD45816_I2`).  
  Preprint: [https://doi.org/10.1101/2024.11.28.625842](https://doi.org/10.1101/2024.11.28.625842)

* **Gene dependency scores** were retrieved from [DepMap](https://depmap.org/portal/data_page/?tab=allData).

* Due to size constraints, full flux sampling data and models are not uploaded to GitHub. They can be generated locally using the notebook `Module1_sampling/script/sampling_script.ipynb`.

* All data required to reproduce the Module 2 analyses of this pipeline are available on Zenodo (DOI: https://doi.org/10.5281/zenodo.17684425).


## Main notebooks description

| Module  | Step | Notebook                                             | Description                                                                                                  |
|--------|------|------------------------------------------------------|--------------------------------------------------------------------------------------------------------------|
| Module 1 | 1️⃣  | sampling_script.ipynb                               | Reconstructs cell specific metabolic models from gene expression data and performs corner based sampling (CBS) for each model. |
| Module 1 | 2️⃣  | comp_nctca_metrics.ipynb                            | Computes the two activity metrics from sampled fluxes. Cycle Propensity and Cycle Flux Intensity.           |
| Module 2 | 1️⃣  | validation_analysis.ipynb                           | Analyzes validation fluxes in ccRCC data and computes proliferation scores for downstream DepMap analysis.  |
| Module 2 | 2️⃣  | data_preprocess.ipynb                               | Preprocesses expression datasets and identifies gene co expression modules.                                 |
| Module 2 | 3️⃣  | elasticnet.ipynb                                    | Performs stable feature selection using multitask Elastic Net regression.                                   |
| Module 2 | 4️⃣  | xgboost.ipynb                                       | Performs stable feature selection using multi output XGBoost models.                                        |
| Module 2 | 5️⃣  | create_subsets.ipynb                                | Builds reduced gene panels from stably selected features.                                                   |
| Module 2 | 6️⃣  | subset_training_shap.ipynb                          | Benchmarks and trains models using reduced gene panels, including SHAP based model interpretation.         |
| Module 2 | 7️⃣  | depmap_analysis.ipynb                               | Correlates gene dependency (DepMap) scores with non canonical TCA cycle activity metrics.                   |
| Module 2 | 8️⃣  | plot.ipynb                                          | Generates the plots and figures used in the manuscript.                                                     |









