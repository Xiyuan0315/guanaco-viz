# Docker Deployment Guide for Guanaco-Viz

## Quick Start

### 1. Build the Image
```bash
docker build -t guanaco-viz .
```

### 2. Run with Docker Compose (Recommended)
```bash
# Update paths in docker-compose.yml first
docker-compose up -d
```

### 3. Access Your App
Visit http://localhost:8080

## Manual Docker Commands

### Basic Run
```bash
docker run -p 8080:8080 guanaco-viz
```

### Run with Data Volume
```bash
docker run -p 8080:8080 \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/config.json:/app/config.json \
  guanaco-viz
```

### Run with Environment Variables
```bash
docker run -p 8080:8080 \
  -e GUANACO_CONFIG=/app/config.json \
  -e GUANACO_DATA_DIR=/app/data \
  -e GUANACO_MAX_CELLS=100000 \
  guanaco-viz
```

## Configuration

### Environment Variables
- `GUANACO_CONFIG`: Path to config.json file (default: /app/config.json)
- `GUANACO_DATA_DIR`: Path to data directory (default: /app/data)
- `GUANACO_MAX_CELLS`: Maximum cells to load (default: 50000)

### Volume Mounts
- `/app/data`: Mount your data directory
- `/app/config.json`: Mount your configuration file

## Production Deployment

### With Docker Compose + Nginx
```bash
docker-compose --profile production up -d
```

### Push to Registry
```bash
# Tag for registry
docker tag guanaco-viz your-registry/guanaco-viz:latest

# Push to registry
docker push your-registry/guanaco-viz:latest
```

## Cloud Deployment Options

### AWS ECS/Fargate
- Use the built image with ECS task definitions
- Mount EFS for persistent data storage

### Google Cloud Run
```bash
# Build and push to GCR
docker tag guanaco-viz gcr.io/PROJECT-ID/guanaco-viz
docker push gcr.io/PROJECT-ID/guanaco-viz

# Deploy to Cloud Run
gcloud run deploy guanaco-viz \
  --image gcr.io/PROJECT-ID/guanaco-viz \
  --platform managed \
  --port 8080
```

### Azure Container Instances
```bash
az container create \
  --resource-group myResourceGroup \
  --name guanaco-viz \
  --image your-registry/guanaco-viz:latest \
  --ports 8080
```

## Troubleshooting

### Container Logs
```bash
docker logs guanaco-viz
```

### Interactive Shell
```bash
docker run -it guanaco-viz /bin/bash
```

### Health Check
```bash
curl http://localhost:8080/
```