from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.decomposition import PCA
import numpy as np

class FixedModuleAggregator(BaseEstimator, TransformerMixin):
    """Aggregate genes into one feature per module by ALWAYS keeping PC1.
    PC1 is fit on TRAIN only within the pipeline; TEST is projected with fixed loadings.
    Input X is expected to be z-scored upstream by StandardScaler in the pipeline.
    """
    def __init__(self, modules, feature_names, random_state=42):
        # Do not transform constructor params (sklearn clone rule)
        self.modules = modules
        self.feature_names = feature_names
        self.random_state = random_state

    def fit(self, X, y=None):
        self.genes_ = np.array(self.feature_names, dtype=object)
        self.loadings_ = []      # list of (idx_array, components[1 x n_genes_in_module])
        self.var_exp_ = []       # PC1 variance explained (NaN for singletons)
        self.out_features_ = []  # module names (one feature per module)

        for mname, g_list in self.modules.items():
            idx = [np.where(self.genes_ == g)[0][0] for g in g_list if g in self.genes_]
            if not idx:
                continue
            idx = np.array(idx, dtype=int)
            Xi = X[:, idx]  # TRAIN (already z-scored by previous scaler)

            if len(idx) == 1:
                # single-gene module: passthrough (identity)
                w = np.ones((1, 1), dtype=np.float32)
                self.loadings_.append((idx, w))
                self.var_exp_.append(np.nan)
                self.out_features_.append(mname)
                continue

            pca = PCA(n_components=1, svd_solver='full', random_state=self.random_state)
            pca.fit(Xi)
            comps = pca.components_.astype(np.float32)  # (1, n_genes_in_module)
            var1  = float(pca.explained_variance_ratio_[0])

            # Sign anchoring: align PC1 with the module mean
            pc1_scores = Xi @ comps.T
            mean_signal = Xi.mean(axis=1, keepdims=True)
            corr = np.corrcoef(pc1_scores.ravel(), mean_signal.ravel())[0, 1]
            if np.isfinite(corr) and corr < 0:
                comps *= -1.0

            self.loadings_.append((idx, comps))
            self.var_exp_.append(var1)
            self.out_features_.append(mname)

        self.out_features_ = np.array(self.out_features_, dtype=object)
        return self

    def transform(self, X):
        blocks = []
        for (idx, comps) in self.loadings_:
            Xi = X[:, idx]
            scores = Xi @ comps.T   # (n_samples, 1)
            blocks.append(scores.astype(np.float32))
        if not blocks:
            return np.empty((X.shape[0], 0), dtype=np.float32)
        return np.hstack(blocks)

    def get_feature_names_out(self, input_features=None):
        return self.out_features_