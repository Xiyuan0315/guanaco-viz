# Guanaco-Viz Deployment Guide

This guide covers deploying the guanaco-viz Dash application to various platforms.

## Overview

The guanaco-viz application is a Dash-based web application that visualizes single-cell data and genome browser tracks. It requires configuration files and data directories to function properly.

## Files Created for Deployment

- `Procfile` - Heroku process configuration
- `runtime.txt` - Python version specification
- `guanaco/wsgi.py` - WSGI entry point with environment setup
- `requirements.txt` - Updated with gunicorn for production

## Platform-Specific Deployment

### Heroku Deployment

#### Prerequisites
- Heroku CLI installed
- Git repository initialized
- Heroku account

#### Steps

1. **Prepare your data and configuration:**
   ```bash
   # Create data directory structure
   mkdir -p custom
   
   # Add your configuration file
   # Create guanaco_v2.json with your dataset configuration
   ```

2. **Create Heroku application:**
   ```bash
   heroku create your-app-name
   ```

3. **Configure environment variables (if needed):**
   ```bash
   # Optional: Set custom paths
   heroku config:set GUANACO_CONFIG=guanaco_v2.json
   heroku config:set GUANACO_DATA_DIR=custom
   heroku config:set GUANACO_MAX_CELLS=5000
   ```

4. **Deploy:**
   ```bash
   git add .
   git commit -m "Deploy guanaco-viz"
   git push heroku main
   ```

#### Important Notes for Heroku:
- **Ephemeral filesystem**: Heroku's filesystem is ephemeral, so uploaded files will be lost on dyno restart
- **File size limits**: Large data files (>500MB) cannot be stored directly in the repo
- **Memory limits**: Consider reducing `GUANACO_MAX_CELLS` for memory-constrained dynos

### Docker Deployment

Create a `Dockerfile`:

```dockerfile
FROM python:3.11.9-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY . .

# Set environment variables
ENV GUANACO_CONFIG=guanaco_v2.json
ENV GUANACO_DATA_DIR=custom
ENV GUANACO_MAX_CELLS=8000

# Expose port
EXPOSE 8000

# Run the application
CMD ["gunicorn", "guanaco.wsgi:application", "--bind", "0.0.0.0:8000", "--workers", "1", "--threads", "2", "--timeout", "300"]
```

Build and run:
```bash
docker build -t guanaco-viz .
docker run -p 8000:8000 -v $(pwd)/data:/app/custom guanaco-viz
```

### AWS Elastic Beanstalk

1. **Create `.ebextensions/python.config`:**
   ```yaml
   option_settings:
     aws:elasticbeanstalk:container:python:
       WSGIPath: guanaco.wsgi:application
     aws:elasticbeanstalk:application:environment:
       GUANACO_CONFIG: guanaco_v2.json
       GUANACO_DATA_DIR: custom
       GUANACO_MAX_CELLS: 8000
   ```

2. **Deploy:**
   ```bash
   eb init
   eb create production
   eb deploy
   ```

## Data Management Strategies

### Small Datasets (< 100MB)
- Include data files directly in the repository
- Place in `custom/` directory
- Reference in `guanaco_v2.json` configuration

### Large Datasets (> 100MB)
- Use cloud storage (S3, Google Cloud Storage)
- Mount volumes in containerized deployments
- Use pre-deployment data synchronization scripts

### Example Configuration for Remote Data

For datasets stored in S3 or similar:

```json
{
  "dataset1": {
    "title": "My Dataset",
    "description": "Example dataset from S3",
    "bucket_urls": [
      "https://my-bucket.s3.amazonaws.com/tracks"
    ],
    "ATAC_name": ["ATAC-seq"],
    "max_height": [100],
    "genome": "hg38"
  }
}
```

## Environment Variables

The application supports these environment variables:

- `GUANACO_CONFIG`: Path to configuration JSON file (default: `guanaco_v2.json`)
- `GUANACO_DATA_DIR`: Directory containing data files (default: `custom`)
- `GUANACO_MAX_CELLS`: Maximum cells to load per dataset (default: `8000`)
- `GUANACO_SEED`: Random seed for cell subsampling (optional)
- `PORT`: Port to run the server on (set automatically by most platforms)

## Configuration File Format

Create a `guanaco_v2.json` file with your dataset configuration:

```json
{
  "color": ["#E69F00", "#56B4E9", "#009E73", "#F0E442"],
  "genome": "hg38",
  "dataset1": {
    "title": "Example Dataset",
    "description": "Description of your dataset",
    "anndata": ["data/dataset1.h5ad"],
    "markers": ["GENE1", "GENE2", "GENE3"],
    "bucket_urls": ["https://bucket.s3.region.amazonaws.com"],
    "ATAC_name": ["ATAC-seq"],
    "max_height": [100]
  }
}
```

## Troubleshooting

### Common Issues

1. **Missing configuration file:**
   - Ensure `guanaco_v2.json` exists in the root directory
   - Check `GUANACO_CONFIG` environment variable

2. **Data files not found:**
   - Verify data directory exists and contains expected files
   - Check `GUANACO_DATA_DIR` environment variable
   - Ensure file paths in configuration are relative to data directory

3. **Memory issues:**
   - Reduce `GUANACO_MAX_CELLS` value
   - Consider using smaller datasets for deployment
   - Monitor application memory usage

4. **Port binding issues:**
   - Ensure the application binds to `0.0.0.0:$PORT`
   - Check that `PORT` environment variable is available

### Logs and Debugging

For Heroku:
```bash
heroku logs --tail
```

For Docker:
```bash
docker logs container-name
```

### Performance Optimization

1. **Reduce cell count**: Lower `GUANACO_MAX_CELLS` for faster loading
2. **Optimize data formats**: Use compressed HDF5 files
3. **Caching**: Consider implementing Redis for session data
4. **CDN**: Use CDN for static assets in production

## Security Considerations

1. **Environment variables**: Never commit sensitive data to version control
2. **Data access**: Ensure S3 buckets have appropriate permissions
3. **HTTPS**: Use HTTPS in production deployments
4. **Authentication**: Consider adding authentication for sensitive datasets

## Monitoring and Maintenance

1. **Health checks**: Implement application health endpoints
2. **Logging**: Configure structured logging for production
3. **Backups**: Regular backups of configuration and small datasets
4. **Updates**: Keep dependencies updated for security

## Support

For deployment issues:
1. Check application logs for detailed error messages
2. Verify all required files and environment variables are present
3. Test locally with production-like configuration
4. Consider the specific constraints of your deployment platform