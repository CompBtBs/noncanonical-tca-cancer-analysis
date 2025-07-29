# noncanonical-tca-cancer-analysis


## Abstract
Cancer cells undergo extensive metabolic rewiring to support growth, survival, and phenotypic plas-
ticity. A non-canonical variant of the tricarboxylic acid (TCA) cycle, characterized by mitochondrial-to-
cytosolic citrate export, has emerged as critical for embryonic stem cell differentiation. However, its role
in cancer remains poorly understood.
Here, we present a two-step computational framework to systematically analyze the activity of this
non-canonical TCA cycle across over 500 cancer cell lines and investigate its role in shaping hallmarks
of malignancy. First, we applied constraint-based modeling to infer the cycle activity, defining two
complementary metrics: Cycle Propensity, reflecting the likelihood of its engagement in each cell line,
and Cycle Flux Intensity, quantifying flux through the pathway’s rate-limiting step. We identified distinct
lineage-specific patterns of pathway utilization. Notably, cells with high Cycle Propensity preferentially
rerouted cytosolic citrate via aconitase 1 (ACO1) and isocitrate dehydrogenase 1 (IDH1), promoting
α-ketoglutarate (αKG) and NADPH production. Elevated engagement of this non-canonical pathway
strongly correlated with Warburg-like metabolic shifts, including decreased oxygen consumption and
increased lactate secretion.
In the second step, to uncover transcriptional correlates of pathway activity, we performed machine
learning–based feature selection using ElasticNet and Random Forest, identifying robust gene signatures
predictive of the non-canonical TCA cycle activity. Over-representation analysis revealed enrichment
in genes involved in invasiveness, angiogenesis, stemness, and key oncogenic pathways. Analysis of
DepMap gene dependency data revealed that non-canonical TCA cycle activity correlates with differ-
ential vulnerability to perturbation of these oncogenic pathways, reinforcing the functional relevance
of identified transcriptional signatures. To further interpret the predictive models, SHapley Additive
exPlanations (SHAP) was applied to prioritize genes contributing most to non-canonical cicle activity,
suggesting novel candidates for experimental investigation.
Overall, our framework enables comprehensive analysis of non-canonical TCA cycle dynamics and
uncovers potential links between metabolic plasticity and malignant phenotypes.
