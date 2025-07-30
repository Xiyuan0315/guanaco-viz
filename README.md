# GUANACO

GUANACO: Interactive visualization tool for single-cell data and genome browser

## Project Structure

```
guanaco-viz/
├── guanaco/                    # Main package directory
│   ├── app.py                  # Dash application initialization
│   ├── cli.py                  # Command-line interface
│   ├── config.py               # Configuration handling
│   ├── data_loader.py          # Data loading utilities
│   ├── layout.py               # Main application layout
│   ├── main.py                 # Application entry point
│   ├── memory_utils.py         # Memory management utilities
│   ├── progress_utils.py       # Progress tracking utilities
│   ├── assets/                 # Static assets (CSS, images)
│   │   ├── scientific_style.css
│   │   ├── logo.png
│   │   └── ...
│   └── pages/                  # Application pages
│       ├── browser/            # Genome browser functionality
│       │   ├── gene_browser.py
│       │   └── utils.py
│       └── single_cell/        # Single-cell analysis pages
│           ├── cellplotly/     # Plotly-based visualizations
│           │   ├── embedding.py
│           │   ├── heatmap1.py
│           │   ├── heatmap2.py
│           │   ├── violin1.py
│           │   ├── violin2.py
│           │   ├── dotmatrix_optimized.py
│           │   ├── pseudotime.py
│           │   └── stacked_bar.py
│           ├── mod021_heatmap.py
│           ├── mod022_violin.py
│           ├── mod023_dotplot.py
│           ├── mod024_stacked_bar.py
│           ├── mod025_pseudotime.py
│           └── single_cell_plots.py
├── example_config.json         # Example configuration file
├── setup.py                    # Package setup
├── requirements.txt            # Python dependencies
└── pixi.toml                   # Pixi package manager config
```

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
- `--seed`: Random seed for cell subsampling

## Configuration

Create a configuration JSON file specifying your datasets. The configuration file should contain:

- **Dataset definitions**: Each dataset with metadata, file paths, and visualization parameters
- **Color schemes**: Custom color palettes for visualizations
- **ATAC-seq tracks**: Genome browser track configurations
- **Marker genes**: Gene lists for specific analyses

See `example_config.json` for a complete example configuration.

## Features

### Single-Cell Analysis
- **Interactive embeddings**: UMAP/t-SNE visualizations with cell metadata overlay
- **Gene expression heatmaps**: Clustered heatmaps with customizable parameters
- **Violin plots**: Distribution plots for gene expression across cell types
- **Dot plots**: Gene expression intensity and percentage visualization
- **Stacked bar charts**: Cell type composition analysis
- **Pseudotime analysis**: Trajectory analysis visualization

### Genome Browser
- **ATAC-seq visualization**: Peak and track display
- **Multi-track support**: Multiple genomic tracks in a single view
- **Interactive navigation**: Zoom and pan functionality
- **Gene annotation**: Integrated gene model display

## License

MIT License
