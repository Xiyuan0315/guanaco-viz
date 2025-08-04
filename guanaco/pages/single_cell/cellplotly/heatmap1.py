import plotly.graph_objects as go
from plotly.subplots import make_subplots
from dash.exceptions import PreventUpdate
import numpy as np
import pandas as pd
from guanaco.data_loader import color_config
from guanaco.pages.single_cell.cellplotly.gene_extraction_utils import (
    extract_gene_expression, extract_multiple_genes, apply_transformation, bin_cells_for_heatmap
)


def plot_heatmap1(adata, genes, labels, adata_obs, groupby, transformation=None, boundary=False, color_map='Viridis', groupby_label_color_map=None, max_cells=50000, n_bins=10000):
    num_genes = len(genes)
    num_labels = len(labels)
    if num_genes == 0 or num_labels == 0:
        raise PreventUpdate

    # Initialize variables for tracking filtered cells
    is_backed = hasattr(adata, 'isbacked') and adata.isbacked
    filtered_obs = adata.obs
    filtered_obs_names = adata.obs_names
    
    # Filter data based on selected labels
    original_adata = adata  # Keep reference to original data
    filtered_adata = None  # Initialize for non-backed case
    
    if labels:
        cell_indices = adata.obs[groupby].isin(labels)
        if is_backed:
            # For backed AnnData, work with indices directly
            cell_indices_array = np.where(cell_indices)[0]
            filtered_obs = adata.obs.iloc[cell_indices_array]
            filtered_obs_names = adata.obs_names[cell_indices_array]
        else:
            # For regular AnnData, create a view
            filtered_adata = adata[cell_indices]
            filtered_obs = filtered_adata.obs
            filtered_obs_names = filtered_adata.obs_names
    else:
        # No filtering - use original data
        filtered_obs = adata.obs
        filtered_obs_names = adata.obs_names
    
    # Check if we need binning for large datasets
    use_binning = len(filtered_obs) > max_cells
    
    genes = list(reversed(genes))
    # Filter out genes that don't exist in the dataset
    valid_genes = [gene for gene in genes if gene in original_adata.var_names]
    if not valid_genes:
        # Return empty figure if no valid genes
        fig = go.Figure()
        fig.add_annotation(
            text="No valid genes found in the dataset",
            xref="paper", yref="paper",
            x=0.5, y=0.5,
            showarrow=False,
            font=dict(size=14)
        )
        fig.update_layout(
            plot_bgcolor='white',
            paper_bgcolor='white',
            height=400
        )
        return fig
    
    # Use centralized gene extraction
    if labels and not is_backed and filtered_adata is not None:
        # For non-backed data, use the filtered adata view we created
        gene_df = extract_multiple_genes(filtered_adata, valid_genes)
    elif labels and is_backed:
        # For backed data with filtering, extract genes individually and filter
        gene_df_list = []
        for gene in valid_genes:
            gene_expr = extract_gene_expression(original_adata, gene)  # Use original adata
            gene_expr_filtered = gene_expr[cell_indices_array] if 'cell_indices_array' in locals() else gene_expr
            gene_df_list.append(pd.Series(gene_expr_filtered, name=gene, index=filtered_obs_names))
        gene_df = pd.concat(gene_df_list, axis=1)
    else:
        # No filtering or fallback - extract all genes from original data
        gene_df = extract_multiple_genes(original_adata, valid_genes)
        
    # Apply transformations using centralized utility
    if transformation:
        for gene in valid_genes:
            gene_df[gene] = apply_transformation(gene_df[gene], method=transformation)
    
    # Add groupby column
    df = gene_df.copy()
    df[groupby] = filtered_obs[groupby].values
    heatmap_df = df.copy()

    # Sort by groupby maintaining the order of labels if provided
    if labels:
        heatmap_df[groupby] = pd.Categorical(
            heatmap_df[groupby],
            categories=labels,
            ordered=True
        )
    sorted_heatmap_df = heatmap_df.sort_values(groupby)

    # Apply binning if needed for large datasets
    if use_binning:
        binned_df = bin_cells_for_heatmap(sorted_heatmap_df, valid_genes, groupby, n_bins)
        sorted_heatmap_df = binned_df

    # Prepare label list and count of cells per group BEFORE creating matrix
    label_list = sorted_heatmap_df[groupby].unique().tolist()
    value_list = [sorted_heatmap_df[sorted_heatmap_df[groupby] == item].shape[0] for item in label_list]
    
    # Create gene matrix with memory optimization
    heatmap_gene_matrix = sorted_heatmap_df[valid_genes].values.T
    if heatmap_gene_matrix.size == 0:
        raise PreventUpdate
    
    del sorted_heatmap_df

    # Use provided color map or create one based on all unique labels from adata_obs
    if groupby_label_color_map is None:
        # Get all unique labels from the original data and sort them for consistent ordering
        unique_labels = sorted(adata_obs[groupby].unique())
        
        # Create consistent color mapping based on all possible labels
        color_map1 = {
            str(label): color_config[i % len(color_config)] for i, label in enumerate(unique_labels)
        }
    else:
        color_map1 = groupby_label_color_map

    y_position = list(range(len(valid_genes)))
    heatmap_height = 40 * len(valid_genes)
    bar_chart_height = 40
    total_height = [heatmap_height, bar_chart_height]
    total_x_range = sum(value_list)

    # Hover text removed for performance - using hoverinfo='skip' instead

    fig = make_subplots(
        rows=2, cols=1,
        row_heights=total_height,
        shared_xaxes=True,
        vertical_spacing=0.02
    )

    # Set colorbar range based on transformation type
    if transformation == "z_score":
        # For z-score, use symmetric range centered at 0
        z_max = max(abs(heatmap_gene_matrix.min()), abs(heatmap_gene_matrix.max()))
        zmin, zmax = -z_max, z_max
        zmid = 0
    else:
        # For other transformations, use data range
        zmin = heatmap_gene_matrix.min()
        zmax = heatmap_gene_matrix.max() 
        zmid = None
    
    colorbar_config = dict(
        title=f'Expression({transformation})',
        len=0.5,
        yanchor='top',
        y=1.0
    )
    
    
    fig.add_trace(go.Heatmap(
        z=heatmap_gene_matrix,
        x=list(range(heatmap_gene_matrix.shape[1])),  # Use matrix shape instead of deleted df
        y=valid_genes,
        colorscale=color_map,
        colorbar=colorbar_config,
        hoverinfo='skip',  # Disable hover
        zmin=zmin,
        zmax=zmax,
        zmid=zmid,
    ), row=1, col=1)
    
    del heatmap_gene_matrix

    if boundary is not False:
        boundary_width = boundary
        boundary_positions = np.cumsum(value_list[:-1]).tolist()
        for pos in boundary_positions:
            fig.add_shape(
                type="line",
                x0=pos, y0=-0.5, x1=pos, y1=len(valid_genes) - 0.5,
                xref='x1', yref='y1',
                line=dict(color="rgba(0, 0, 0, 0.8)", width=boundary_width),
                row=1, col=1
            )

    x_pos = 0
    for i, label in enumerate(label_list):
        fig.add_trace(go.Bar(
            x=[value_list[i]],
            marker_color=color_map1[str(label)],
            name=str(label),
            hoverinfo='skip',  # Disable hover
            showlegend=False,
            orientation='h',
            width=1
        ), row=2, col=1)

        fig.add_annotation(
            x=x_pos + value_list[i] / 2,
            y=-0.5,
            text=str(label),
            showarrow=False,
            xanchor='center',
            yanchor='top',
            font=dict(size=10, color='black'),
            textangle=-90,
            row=2, col=1
        )
        x_pos += value_list[i]

    fig.update_layout(
        barmode='stack',
        showlegend=False,
        plot_bgcolor='white',
        paper_bgcolor='white',
        xaxis=dict(
            range=[0, total_x_range],
            constrain='domain',
            constraintoward='center'
        ),
        xaxis2=dict(
            range=[0, total_x_range], 
            visible=False,
            constrain='domain',
            constraintoward='center'
        ),
        yaxis=dict(
            tickmode='array',
            tickvals=y_position,
            ticktext=valid_genes,
            tickfont=dict(size=14),
            constrain='domain',
            constraintoward='middle'
        ),
        yaxis2=dict(
            tickmode='array',
            tickvals=[0],
            ticktext=[groupby],
            showgrid=False,
            tickfont=dict(size=14, family='Arial Black', color='black'),
            constrain='domain',
            constraintoward='middle'
        ),
        height=max(450, sum(total_height)),
        margin=dict(t=50, b=150, l=50, r=50),
        hovermode='closest'
    )
    fig.update_xaxes(tickangle=90, visible=False)

    return fig
