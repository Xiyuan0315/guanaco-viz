# GUANACO

GUANACO: Interactive visualization tool for single-cell data and genome browser

## Installation

```bash
pip install guanaco-viz
```

## Usage

```bash
guanaco -c config.json
```

### Command-line Options

- `-c, --config`: Path to configuration JSON file (default: guanaco_v2.json)
- `-d, --data-dir`: Directory containing AnnData files (default: custom)
- `-p, --port`: Port to run the server on (default: 4399)
- `--host`: Host to run the server on (default: 0.0.0.0)
- `--debug`: Run server in debug mode
- `--max-cells`: Maximum number of cells to load per dataset (default: 8000)
- `--seed`: Random seed for cell subsampling

## Configuration

Create a configuration JSON file specifying your datasets. See `example_config.json` for an example.

## License

MIT License