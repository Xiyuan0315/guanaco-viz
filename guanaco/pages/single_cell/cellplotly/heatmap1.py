import plotly.graph_objects as go
from plotly.subplots import make_subplots
from dash.exceptions import PreventUpdate
import numpy as np
import pandas as pd
from guanaco.data_loader import color_config

def plot_heatmap1(adata, genes, labels, adata_obs, groupby, transformation=None, boundary=False, color_map='Viridis', groupby_label_color_map=None):
    num_genes = len(genes)
    num_labels = len(labels)
    if num_genes == 0 or num_labels == 0:
        raise PreventUpdate

    # Filter data based on selected labels
    if labels:
        cell_indices = adata.obs[groupby].isin(labels)
        adata = adata[cell_indices]
    
    genes = list(reversed(genes))
    # Filter out genes that don't exist in the dataset
    valid_genes = [gene for gene in genes if gene in adata.var_names]
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
    
    # For backed AnnData, we need to extract data directly without creating views
    if hasattr(adata, 'isbacked') and adata.isbacked:
        # Extract gene indices for valid genes
        gene_indices = [adata.var_names.get_loc(gene) for gene in valid_genes]
        # Extract expression data directly from the backed file
        if hasattr(adata.X, 'toarray'):
            # For sparse backed data
            gene_expression_matrix = adata.X[:, gene_indices].toarray()
        else:
            # For dense backed data
            gene_expression_matrix = adata.X[:, gene_indices]
    else:
        # Original code for non-backed AnnData
        adata_selected = adata[:, valid_genes]
        if hasattr(adata_selected.X, 'toarray'):
            gene_expression_matrix = adata_selected.X.toarray()
        else:
            gene_expression_matrix = adata_selected.X

    if transformation == "log":
        gene_expression_matrix = np.log1p(gene_expression_matrix)
    elif transformation == "z_score":
        gene_expression_matrix = (gene_expression_matrix - gene_expression_matrix.mean(axis=0)) / gene_expression_matrix.std(axis=0)

    df = pd.DataFrame(gene_expression_matrix, columns=valid_genes, index=adata.obs_names)
    df[groupby] = adata.obs[groupby].values
    heatmap_df = df.copy()

    # Sort by groupby maintaining the order of labels if provided
    if labels:
        heatmap_df[groupby] = pd.Categorical(
            heatmap_df[groupby],
            categories=labels,
            ordered=True
        )
    sorted_heatmap_df = heatmap_df.sort_values(groupby)

    heatmap_gene_matrix = sorted_heatmap_df[valid_genes].values.T
    if heatmap_gene_matrix.size == 0:
        raise PreventUpdate

    # Prepare label list and count of cells per group
    label_list = sorted_heatmap_df[groupby].unique().tolist()
    value_list = [sorted_heatmap_df[sorted_heatmap_df[groupby] == item].shape[0] for item in label_list]

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

    # Hover text for the heatmap
    hover_text = np.array([
        [f"Gene: {gene}<br>{groupby}: {group}<br>Expression: {expr:.2f}"
         for group, expr in zip(sorted_heatmap_df[groupby], row)]
        for gene, row in zip(valid_genes, heatmap_gene_matrix)
    ])

    fig = make_subplots(
        rows=2, cols=1,
        row_heights=total_height,
        shared_xaxes=True,
        vertical_spacing=0.02
    )

    fig.add_trace(go.Heatmap(
        z=heatmap_gene_matrix,
        x=list(range(len(sorted_heatmap_df[groupby]))),
        y=valid_genes,
        colorscale=color_map,
        colorbar=dict(
            title=f'Expression({transformation})',
            len=0.5,
            yanchor='top',
            y=1.0
        ),
        text=hover_text,
        hoverinfo='text',
        zmin=heatmap_gene_matrix.min(),
        zmax=heatmap_gene_matrix.max(),
    ), row=1, col=1)

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
            hovertext=f'{value_list[i]} ({label})',
            hoverinfo='text',
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
