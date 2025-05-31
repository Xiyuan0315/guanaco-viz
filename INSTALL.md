# GUANACO Installation Guide

## Installation from PyPI (when published)

```bash
pip install guanaco-viz
```

## Installation from Source

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

### 3. Install from distribution files
If you have the wheel file:
```bash
pip install dist/guanaco_viz-0.1.0-py3-none-any.whl
```

## Running GUANACO

After installation, you can run GUANACO using:

```bash
guanaco -c your_config.json
```

See `guanaco --help` for all available options.

## Publishing to PyPI

To publish to PyPI (for maintainers):

1. Install build and twine:
```bash
pip install build twine
```

2. Build the distribution:
```bash
python -m build
```

3. Upload to PyPI:
```bash
python -m twine upload dist/*
```

## Dependencies

All dependencies will be automatically installed. See `requirements.txt` for the full list.