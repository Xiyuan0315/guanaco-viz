# GUANACO
GUANACO (Graphical Unified Analysis and Navigation of Cellular Omics), a Python-based platform tailored for biologists to explore multi-omics single-cell data.<img width="100" height="140" alt="logo" src="https://github.com/user-attachments/assets/d605e945-7fa1-493c-8444-536b55c15ff4" />

<img width="2009" height="2252" alt="Copy of Copy of Copy of figure1_v3 (1)" src="https://github.com/user-attachments/assets/3f3b5557-4605-452e-9aea-97cda45e35fd" />



## Installation

### 1. Clone the repository
```bash
git clone https://github.com/Systems-Immunometabolism-Lab/GUANACO_updated
cd GUANACO_updated
```

### 2. Install from local directory
```bash
pip install .
```
Or for development (editable install):
```bash
pip install -e .
```
## Usage
```bash
guanaco -c config.json -d data_folder
```

### Command-line Options

- `-c, --config`: Name of configuration JSON file (relative to --data-dir) (default: guanaco.json)
- `-d, --data-dir`: Directory containing AnnData files referenced in config (default: current directory)
- `-p, --port`: Port to run the Dash server on (default: 4399)
- `--host`: Host to run the Dash server on (default: 0.0.0.0)
- `--debug`: Run server in debug mode
- `--max-cells`: Maximum number of cells to load per dataset (default: 8000)
- `--backed-mode`: Enable backed mode for memory-efficient loading of large datasets

## Configuration

Create a configuration JSON file specifying your datasets. See `example_config.json` for a complete example configuration. Simpliest case for visualizing scRNA data(.h5ad) is:
```
{
  "Demo": {"sc_data": "PBMC_int.h5ad"}
}
```
## Features

### Single-Cell Analysis
- **Interactive embeddings**: Linked UMAP/t-SNE visualizations with cell metadata/gene expression overlay. Select cells from embeddings to 
- **Gene expression heatmaps**: Clustered heatmaps with customizable parameters
- **Violin plots**: Distribution plots for gene expression across cell types
- **Dot plots**: Gene expression intensity and percentage visualization
- **Stacked bar charts**: Cell type composition analysis
- **Pseudotime analysis**: Trajectory analysis visualization

### Genome Browser
- **ATAC-seq visualization**: Peak and track display
- **Multi-track support**: Multiple genomic tracks in a single view
- **Interactive navigation**: Zoom and pan
- **Motif Database**: Jaspar database for searching

