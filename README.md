# GPT Cell Annotator

This package provides a wrapper around [scanpy](https://github.com/theislab/scanpy) and [OpenAI GPT](https://platform.openai.com/docs/) to perform cell-type annotation on single-cell RNA-seq data.

## Installation

```bash
pip install -e "git+https://github.com/VPetukhov/GPTCellAnnotator.git"
```

## Usage

### Set OpenAI key

The package requires you to have an OpenAI API key. You can get one [here](https://help.openai.com/en/articles/4936850-where-do-i-find-my-secret-api-key). Once you have it, you can set it in your notebook to `openai.api_key`. You can save it to a `.env` file of your project and then use [`python-dotenv`](https://pypi.org/project/python-dotenv/) to read it to an environment variable:

```python
from dotenv import load_dotenv
load_dotenv()

openai.api_key = os.getenv("OPENAI_API_KEY")
```

### Package workflow

1. **Generate a list of expected cell types and their markers** for the given tissue using GPT. This list would later be used in the prompt to improve the quality. First, it is needed to make cell type names standardized. Second, it makes GPT focus on the cell types that are relevant to the tissue of interest.
2. **Generate a list of markers for the given cell types** using GPT. GPT is surprisingly good in providing relevant cell type markers. While this list would not be directly used for cell type annotation, it seems to improve the annotation quality if provided to the annotation prompt. These two steps can be run with:
    ```python
    expected_types, expected_markers = get_expected_cell_types(
        species='mouse', tissue='pancreas', model='gpt-4', max_tokens=800
    )
    ```
3. **Process scRNA-seq data** using `scanpy`. This step includes filtering, normalization, dimensionality reduction, clustering, and marker gene selection. In addition to the standard `scanpy` pipeline, it uses `AUC` and `Specificity` metrics to rank genes. Example processing for the `AnnData` object `ad`:
    ```python
    sc.pp.neighbors(ad, n_neighbors=30)
    sc.tl.leiden(ad, resolution=0.8)
    marker_dfs = get_markers_per_cluster(ad, clustering='leiden')
    ```
4. **Annotate the scRNA-seq clusters based on their markers** using GPT. This step takes the data from the three previous steps, iterates over the clusters and queries OpenAI for each of them.
    ```python
    annotation_res = annotate_clusters(
        marker_genes, species='mouse', tissue='pancreas', expected_markers=expected_markers,
        model='gpt-4'
    )
    ann_df = parse_annotation(annotation_res)
    ```

After these steps, the data can be visualized with:
```python
ad.obs['annotation'] = ann_df['Cell type'][ad.obs['leiden'].values.astype(int)].values
sc.pl.umap(ad, color=['annotation'], legend_loc='on data', legend_fontsize=8)
```


## Data privacy

Biological data is usually sensitive and should be handled with care. This package does not send any actual data to OpenAI. All processing goes locally, and the only data that is sent to OpenAI are:
1. Name of the tissue and the species
2. List of clusters and top-5 marker genes per cluster (only names, not expression values)
