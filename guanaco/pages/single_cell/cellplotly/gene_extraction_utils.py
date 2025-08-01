import numpy as np
import pandas as pd
from scipy.sparse import issparse
from functools import lru_cache
import hashlib
import time
import gc


class GeneExpressionCache:
    """Cache for gene expression data to avoid repeated extractions."""
    
    def __init__(self, max_size=64, max_age_seconds=1800):  # Reduced to 64 items, 30-min expiration
        self.cache = {}
        self.max_size = max_size
        self.max_age_seconds = max_age_seconds
        self.timestamps = {}
        self.access_order = []
        
    def _make_key(self, adata, gene, layer=None, use_raw=False):
        """Create a unique key for caching."""
        # For backed mode, use the filename as a stable identifier
        if hasattr(adata, 'isbacked') and adata.isbacked and hasattr(adata, 'filename'):
            adata_id = adata.filename
        else:
            adata_id = id(adata)
        
        # Include data shape to differentiate between filtered and unfiltered data
        data_shape = adata.shape
        
        return (adata_id, gene, layer, use_raw, data_shape)
    
    def get(self, adata, gene, layer=None, use_raw=False):
        """Get gene expression from cache or compute if not cached."""
        key = self._make_key(adata, gene, layer, use_raw)
        current_time = time.time()
        
        # Check if cached and not expired
        if key in self.cache:
            if current_time - self.timestamps.get(key, 0) < self.max_age_seconds:
                # Move to end of access order
                if key in self.access_order:
                    self.access_order.remove(key)
                self.access_order.append(key)
                return self.cache[key]  # Return reference, not copy (saves memory)
            else:
                # Expired - remove from cache
                del self.cache[key]
                del self.timestamps[key]
                if key in self.access_order:
                    self.access_order.remove(key)
        
        # Compute if not in cache - use extract without cache to avoid recursion
        expr = extract_gene_expression(adata, gene, layer, use_raw, use_cache=False)
        
        # Add to cache with size limit
        if len(self.cache) >= self.max_size:
            # Remove least recently used
            if self.access_order:
                lru_key = self.access_order.pop(0)
                del self.cache[lru_key]
                del self.timestamps[lru_key]
        
        self.cache[key] = expr
        self.timestamps[key] = current_time
        self.access_order.append(key)
        
        # Periodic garbage collection
        if len(self.cache) % 20 == 0:
            gc.collect()
        
        return expr
    
    def clear(self):
        """Clear the cache and free memory."""
        self.cache.clear()
        self.timestamps.clear()
        self.access_order.clear()
        gc.collect()


def extract_gene_expression(adata, gene, layer=None, use_raw=False, use_cache=True):
    """
    Optimized gene expression extraction with caching support.
    
    Parameters:
    -----------
    adata : AnnData
        The annotated data object
    gene : str
        Gene name to extract
    layer : str, optional
        Layer to extract from (default: None uses .X)
    use_raw : bool
        Whether to use raw data (default: False)
    use_cache : bool
        Whether to use caching (default: True)
    
    Returns:
    --------
    np.ndarray
        Gene expression values as 1D array
    """
    # Use cache if enabled
    if use_cache:
        return _gene_cache.get(adata, gene, layer, use_raw)
    
    # Direct extraction without cache
    # Use raw data if requested
    data_source = adata.raw if use_raw and adata.raw is not None else adata
    
    # Get gene index
    if gene not in data_source.var_names:
        raise ValueError(f"Gene '{gene}' not found in data")
    
    gene_idx = data_source.var_names.get_loc(gene)
    
    # Get expression matrix
    if layer is not None:
        if layer not in data_source.layers:
            raise ValueError(f"Layer '{layer}' not found")
        X = data_source.layers[layer]
    else:
        X = data_source.X
    
    # Extract gene expression
    gene_data = X[:, gene_idx]
    
    # Handle sparse matrices efficiently
    if issparse(gene_data):
        # Only convert the specific column to dense
        gene_expr = gene_data.toarray().flatten()
    else:
        gene_expr = gene_data.flatten() if hasattr(gene_data, 'flatten') else np.asarray(gene_data).flatten()
    
    return gene_expr


def extract_multiple_genes(adata, genes, layer=None, use_raw=False):
    """
    Extract multiple genes efficiently in a single operation.
    
    Parameters:
    -----------
    adata : AnnData
        The annotated data object
    genes : list of str
        List of gene names to extract
    layer : str, optional
        Layer to extract from (default: None uses .X)
    use_raw : bool
        Whether to use raw data (default: False)
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with genes as columns and cells as rows
    """
    # Use raw data if requested
    data_source = adata.raw if use_raw and adata.raw is not None else adata
    
    # Get gene indices
    gene_indices = []
    valid_genes = []
    for gene in genes:
        if gene in data_source.var_names:
            gene_indices.append(data_source.var_names.get_loc(gene))
            valid_genes.append(gene)
        else:
            print(f"Warning: Gene '{gene}' not found in data")
    
    if not gene_indices:
        raise ValueError("No valid genes found")
    
    # Get expression matrix
    if layer is not None:
        if layer not in data_source.layers:
            raise ValueError(f"Layer '{layer}' not found")
        X = data_source.layers[layer]
    else:
        X = data_source.X
    
    # Extract all genes at once
    gene_data = X[:, gene_indices]
    
    # Convert to dense if sparse
    if issparse(gene_data):
        gene_data = gene_data.toarray()
    elif hasattr(gene_data, 'compute'):  # Handle dask arrays
        gene_data = gene_data.compute()
    
    # Create DataFrame
    df = pd.DataFrame(gene_data, columns=valid_genes, index=adata.obs_names)
    
    return df


def apply_transformation(expr, method='log1p', copy=True):
    """
    Apply transformation to expression data efficiently.
    
    Parameters:
    -----------
    expr : np.ndarray or pd.DataFrame
        Expression data
    method : str
        Transformation method: 'log1p', 'log', 'sqrt', 'zscore', 'z_score', 'scale'
    copy : bool
        Whether to copy data before transformation
    
    Returns:
    --------
    Transformed expression data
    """
    if copy:
        expr = expr.copy()
    
    # Handle both 'zscore' and 'z_score'
    if method in ['zscore', 'z_score']:
        # Z-score normalization
        if isinstance(expr, pd.DataFrame):
            return expr.apply(lambda x: (x - x.mean()) / x.std(), axis=0)
        else:
            mean = np.mean(expr)
            std = np.std(expr)
            return (expr - mean) / std if std > 0 else expr - mean
    elif method == 'log1p' or method == 'log':
        return np.log1p(expr)
    elif method == 'sqrt':
        return np.sqrt(expr)
    elif method == 'scale':
        # Min-max scaling to [0, 1]
        if isinstance(expr, pd.DataFrame):
            return expr.apply(lambda x: (x - x.min()) / (x.max() - x.min()), axis=0)
        else:
            min_val = np.min(expr)
            max_val = np.max(expr)
            return (expr - min_val) / (max_val - min_val) if max_val > min_val else expr - min_val
    else:
        return expr


def get_expression_percentiles(expr, percentiles=[0, 25, 50, 75, 95, 99, 100]):
    """
    Get percentile values for expression data.
    
    Useful for setting vmin/vmax in plots.
    """
    return np.percentile(expr[expr > 0] if len(expr[expr > 0]) > 0 else expr, percentiles)


# Global cache instance
_gene_cache = GeneExpressionCache()


def clear_gene_cache():
    """Clear the global gene expression cache."""
    _gene_cache.clear()


def get_cache_info():
    """Get information about the current cache state."""
    return {
        'size': len(_gene_cache.cache),
        'max_size': _gene_cache.max_size,
        'keys': list(_gene_cache.cache.keys())
    }