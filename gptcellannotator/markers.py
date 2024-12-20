from sklearn.metrics import roc_auc_score
import numpy as np
import scanpy as sc
from pandas import Series,DataFrame
from tqdm.auto import tqdm


def get_specificity(genes: str, clust_mask: np.array, ad: sc.AnnData):
    expr_mask = (ad.raw[:, genes].X > 0).toarray()
    fpr = np.sum(expr_mask & np.atleast_2d(~clust_mask).T, axis=0) / np.sum(~clust_mask)
    return 1 - fpr


def get_auc(genes: str, clust_mask: np.array, ad: sc.AnnData):
    return np.array([roc_auc_score(clust_mask, ad.raw[:, g].X.toarray()[:,0]) for g in genes])


def get_markers_per_cluster(
        ad: sc.AnnData, clustering: str = 'leiden',
        min_specificity: float = 0.75, min_markers: int = 15, max_markers: int = 200
    ):
    sc.tl.rank_genes_groups(ad, groupby=clustering, method="wilcoxon")

    marker_dfs = {}
    for cli in tqdm(ad.obs[clustering].unique()):
        genes = ad.uns['rank_genes_groups']['names'][cli][:max_markers]
        clust_mask = (ad.obs[clustering] == cli)
        specificity = get_specificity(genes, clust_mask, ad)
        mask = (specificity >= min(min_specificity, sorted(specificity)[-min_markers]))
        genes = genes[mask]
        specificity = specificity[mask]
        auc = get_auc(genes, clust_mask, ad)
        marker_dfs[cli] = DataFrame({'gene': genes, 'specificity': specificity, 'auc': auc})

    return Series(marker_dfs)


def get_top_markers(marker_dfs: Series, min_auc: float = 0.7, max_markers: int = 7):
    marker_genes = marker_dfs.map(
        lambda x: x[x.auc > min_auc].sort_values('auc', ascending=False).head(max_markers).gene.values
    )
    return marker_genes[np.argsort(marker_genes.index.astype(int))]
