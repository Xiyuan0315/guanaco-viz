from .dotmatrix_optimized import plot_dot_matrix_optimized

def plot_dot_matrix(adata, genes, groupby, selected_labels, aggregation='mean', transformation=None,standardization = None, vmin=None, vmax=None, expression_threshold=0, color_map='Viridis', plot_type='dotplot'):
    # Always use optimized implementation for better performance and caching
    return plot_dot_matrix_optimized(adata, genes, groupby, selected_labels, aggregation, transformation, standardization, vmin, vmax, expression_threshold, color_map, plot_type)
