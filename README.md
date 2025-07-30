# Mechanistically Informed Machine Learning Links Non-Canonical TCA Cycle Activity to Warburg Metabolism and Hallmarks of Malignancy

This repository contains the code and relevant resources used in the analysis presented in the paper:

**Mechanistically Informed Machine Learning Links Non-Canonical TCA Cycle Activity to Warburg Metabolism and Hallmarks of Malignancy**  
*Lin Lihao, Francesco Lapi, Bruno G. Galuzzi, Marco Vanoni, Lilia Alberghina, Chiara Damiani*  
(Submitted / Preprint available at: [insert DOI or arXiv link])

## Paper abstract

Cancer cells undergo extensive metabolic rewiring to support growth, survival, and phenotypic plasticity. A non-canonical variant of the tricarboxylic acid (TCA) cycle, characterized by mitochondrial-to-cytosolic citrate export, has emerged as critical for embryonic stem cell differentiation. However, its role in cancer remains poorly understood.

Here, we present a two-step computational framework to systematically analyze the activity of this non-canonical TCA cycle across over 500 cancer cell lines and investigate its role in shaping hallmarks of malignancy. First, we applied constraint-based modeling to infer the cycle activity, defining two complementary metrics: *Cycle Propensity*, reflecting the likelihood of its engagement in each cell line, and *Cycle Flux Intensity*, quantifying flux through the pathway’s rate-limiting step. We identified distinct lineage-specific patterns of pathway utilization. Notably, cells with high Cycle Propensity preferentially rerouted cytosolic citrate via aconitase 1 (ACO1) and isocitrate dehydrogenase 1 (IDH1), promoting α-ketoglutarate (αKG) and NADPH production. Elevated engagement of this non-canonical pathway strongly correlated with Warburg-like metabolic shifts, including decreased oxygen consumption and increased lactate secretion.

In the second step, to uncover transcriptional correlates of pathway activity, we performed machine learning–based feature selection using ElasticNet and Random Forest, identifying robust gene signatures predictive of the non-canonical TCA cycle activity. Over-representation analysis revealed enrichment in genes involved in invasiveness, angiogenesis, stemness, and key oncogenic pathways. Analysis of DepMap gene dependency data revealed that non-canonical TCA cycle activity correlates with differential vulnerability to perturbation of these oncogenic pathways, reinforcing the functional relevance of identified transcriptional signatures. To further interpret the predictive models, SHapley Additive exPlanations (SHAP) was applied to prioritize genes contributing most to non-canonical cicle activity, suggesting novel candidates for experimental investigation.

Overall, our framework enables comprehensive analysis of non-canonical TCA cycle dynamics and uncovers potential links between metabolic plasticity and malignant phenotypes.

## Pipeline overview

<p align="center">
  <img src="Graphical_Abstract.png" alt="Graphical Abstract" width="700"/>
</p>

We adopted a two-step computational strategy to dissect non canonical TCA cycle (Arnold Cycle) activity across cancer cell lines.

In the first module, we used constraint-based metabolic modeling to infer the activity of the Arnold Cycle from gene expression data. Specifically, we extracted transcriptomic profiles from 513 cancer cell lines cultured under standardized conditions and reconstructed a cell-specific metabolic model for each line based on the ENGRO2.2 core network. Reaction constraints were individualized by computing Reaction Activity Scores (RAS) using Gene–Protein–Reaction (GPR) rules, and these scores were used to rescale flux boundaries to reflect the transcriptional context of each cell line. To comprehensively explore the space of feasible metabolic states, we performed unbiased flux sampling on each transcriptome-constrained model using a corner-based sampling algorithm. This allowed us to approximate the distribution of steady-state flux configurations consistent with mass balance and expression constraints, without assuming predefined cellular objectives. From these sampled flux distributions, we derived two complementary metrics to quantify Arnold Cycle activity: `Cycle Propensity`, defined as the fraction of flux samples in which the three reactions (*SLC25A1*, *ACLY*, *MDH1*) of the Arnold Cycle are simultaneously active; and `Cycle Flux Intensity`, which captures the average flux through the bottleneck reaction when the cycle is active.

In the second module, these *in silico*–derived metrics were used as training labels for supervised machine learning models. To predict `Cycle Propensity` and `Cycle Flux Intensity` from transcriptional features beyond the metabolic network, we trained regressors using expression data from non-metabolic genes. After model training, we performed feature selection to identify robust predictors of Arnold Cycle activity, thereby uncovering broader transcriptional programs associated with its engagement.


## Data availability

* **Gene expression data** were obtained from the [Cancer Cell Line Encyclopedia (CCLE)](https://sites.broadinstitute.org/ccle/datasets).
* **Gene dependency scores** were retrieved from [DepMap](https://depmap.org/portal/data_page/?tab=allData).
* **Flux metrics** (Cycle Propensity and Cycle Flux Intensity) for each cell line are available as a CSV file in: `Module2_ml/data/info/non_canonical_state.csv`
* Mean flux values across the sampled steady-state distributions are available in:
`Module1_sampling/dati_sampling/dati/sampling/mean`
*  Due to size constraints, full flux sampling data and models are not uploaded to GitHub. If needed, they can be generated locally using thenotebook `Module1_sampling/script/sampling_script.ipynb`.
*  Preprocessed features and labels used for model training are available in:
`Module2_ml/data/X_y/`

## Main notebooks description

| Step | Notebook                                           | Description                                                                                                                                                    |
| ---- | -------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1️⃣  | `Module1_sampling/script/sampling_script.ipynb`    | Reconstructs **cell-specific metabolic models** from gene expression data and performs **corner-based sampling** (CBS) for each model.                         |
| 2️⃣  | `Module1_sampling/script/comp_nctca_metrics.ipynb` | Computes the two activity metrics from sampled fluxes: `Cycle Propensity` and `Cycle Flux Intensity`.                                                          |
| 3️⃣  | `Module2_ml/scripts/en_rf_selec.ipynb`             | Performs **feature selection** using ElasticNet and Random Forest to identify transcriptional predictors of Arnold Cycle activity.                             |
| 4️⃣  | `Module2_ml/scripts/xgboost_shap.ipynb`            | Trains an **XGBoost model** on the selected features and applies **SHAP analysis** to interpret model predictions and prioritize key transcriptional features. |








