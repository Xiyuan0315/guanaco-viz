import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
import plotly.express as px
from guanaco.pages.single_cell.cellplotly.gene_extraction_utils import (
    extract_gene_expression, extract_multiple_genes, apply_transformation, bin_cells_for_heatmap
)

def is_continuous_annotation(adata, annotation, threshold=50):
    """Check if an annotation is continuous based on unique value count and data type."""
    if annotation not in adata.obs.columns:
        return False
    
    # Check data type
    dtype = adata.obs[annotation].dtype
    if dtype in ['float32', 'float64', 'int32', 'int64']:
        # Numeric type - check unique values
        n_unique = adata.obs[annotation].nunique()
        return n_unique >= threshold
    return False

def plot_unified_heatmap(
    adata, genes, groupby1, groupby2=None, labels=None, log=False, z_score=False, 
    boundary=False, color_map='Viridis', groupby1_label_color_map=None, 
    groupby2_label_color_map=None, max_cells=50000, n_bins=10000, transformation=None
):
    """
    Unified heatmap function that handles both single and dual annotation cases.
    
    Parameters:
    - adata: AnnData object
    - genes: List of genes to plot
    - groupby1: Primary annotation column
    - groupby2: Secondary annotation column (optional)
    - labels: Selected labels for filtering (optional)
    - log: Apply log transformation
    - z_score: Apply z-score transformation
    - boundary: Whether to draw boundary lines
    - color_map: Colormap for the heatmap
    - groupby1_label_color_map: Color map for primary annotation
    - groupby2_label_color_map: Color map for secondary annotation
    - max_cells: Maximum cells for binning
    - n_bins: Number of bins for large datasets
    - transformation: Transformation method (for heatmap1 compatibility)
    """
    
    # Handle parameter compatibility between old heatmap1 and heatmap2
    if transformation:
        if transformation == 'log':
            log = True
        elif transformation in ['z_score', 'zscore']:
            z_score = True
    
    # Check if we need to handle continuous secondary annotation
    if groupby2 and is_continuous_annotation(adata, groupby2):
        return plot_heatmap2_continuous(
            adata, genes, groupby1, groupby2, labels, log, z_score,
            color_map, groupby1_label_color_map, max_cells, n_bins
        )
    
    # Initialize variables for tracking filtered cells
    is_backed = hasattr(adata, 'isbacked') and adata.isbacked
    filtered_obs = adata.obs
    filtered_obs_names = adata.obs_names
    
    # Filter data based on selected labels
    original_adata = adata  # Keep reference to original data
    
    if labels:
        cell_indices = adata.obs[groupby1].isin(labels)
        if is_backed:
            cell_indices_array = np.where(cell_indices)[0]
            filtered_obs = adata.obs.iloc[cell_indices_array]
            filtered_obs_names = adata.obs_names[cell_indices_array]
        else:
            adata = adata[cell_indices].copy()
            filtered_obs = adata.obs
            filtered_obs_names = adata.obs_names

    # Filter out genes that don't exist in the dataset
    valid_genes = [gene for gene in genes if gene in original_adata.var_names]
    if not valid_genes:
        fig = go.Figure()
        fig.add_annotation(
            text="No valid genes found in the dataset",
            xref="paper", yref="paper",
            x=0.5, y=0.5,
            showarrow=False,
            font=dict(size=14)
        )
        fig.update_layout(plot_bgcolor='white', paper_bgcolor='white', height=400)
        return fig

    # Check if we need binning for large datasets
    use_binning = len(filtered_obs) > max_cells
    
    # Get expression data
    if labels and is_backed:
        gene_df_list = []
        for gene in valid_genes:
            gene_expr = extract_gene_expression(original_adata, gene)
            gene_expr_filtered = gene_expr[cell_indices_array]
            gene_df_list.append(pd.Series(gene_expr_filtered, name=gene, index=filtered_obs_names))
        gene_df = pd.concat(gene_df_list, axis=1)
    else:
        gene_df = extract_multiple_genes(adata, valid_genes)

    # Apply transformations
    if log:
        for gene in valid_genes:
            gene_df[gene] = apply_transformation(gene_df[gene], method='log1p')
    if z_score:
        for gene in valid_genes:
            gene_df[gene] = apply_transformation(gene_df[gene], method='z_score')

    # Create combined dataframe
    gene_df.insert(0, 'CellID', filtered_obs_names)
    
    # Prepare annotation columns
    annotation_columns = [groupby1]
    if groupby2 and groupby2 != 'None' and groupby2 != groupby1:
        annotation_columns.append(groupby2)
    
    label_df = pd.DataFrame(filtered_obs[annotation_columns])
    label_df.insert(0, 'CellID', filtered_obs_names)
    heatmap_df = pd.merge(gene_df, label_df, on='CellID')

    # Sort data
    if labels:
        heatmap_df[groupby1] = pd.Categorical(heatmap_df[groupby1], categories=labels, ordered=True)
    
    sort_columns = [groupby1]
    if len(annotation_columns) > 1:
        sort_columns.append(groupby2)
    sorted_heatmap_df = heatmap_df.sort_values(sort_columns)

    # Create gene matrix
    heatmap_gene_matrix = sorted_heatmap_df[valid_genes].values.T

    # Calculate label statistics for primary annotation
    label_list1 = sorted_heatmap_df[groupby1].unique().tolist()
    value_list1 = []
    for item in label_list1:
        value_list1.append(sorted_heatmap_df[sorted_heatmap_df[groupby1] == item].shape[0])

    # Setup color maps and determine layout
    has_secondary = len(annotation_columns) > 1
    
    # Import the consistent color config used by other plots
    from guanaco.data_loader import color_config
    
    # Create default color maps using the same colors as other plots
    tab20_colors = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', 
                    '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', 
                    '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', 
                    '#17becf', '#9edae5']
    
    # Get ALL unique labels from the original dataset (not just filtered ones)
    # This ensures consistent color assignment regardless of filtering
    all_unique_labels_primary = sorted(original_adata.obs[groupby1].unique())
    
    # Use the consistent color_config for primary annotation (same as other plots)
    # Assign colors based on ALL labels to maintain consistency
    default_color_map1 = {
        label: color_config[i % len(color_config)] for i, label in enumerate(all_unique_labels_primary)
    }
    
    # Use provided color maps or defaults
    if groupby1_label_color_map:
        # Use provided color map (it should already contain all labels)
        color_map1 = groupby1_label_color_map
    else:
        # Use default colors when no custom color map provided
        color_map1 = default_color_map1

    # Setup for secondary annotation if present
    if has_secondary:
        sorted_heatmap_df['combined'] = sorted_heatmap_df[groupby1].astype(str) + '_' + sorted_heatmap_df[groupby2].astype(str)
        label_list2 = sorted_heatmap_df['combined'].unique().tolist()
        label2_dict = {item: item.split('_')[-1] for item in label_list2}
        value_list2 = []
        for item in label_list2:
            value_list2.append(sorted_heatmap_df[sorted_heatmap_df['combined'] == item].shape[0])
        
        # Get ALL unique secondary labels from original dataset
        all_unique_labels_secondary = sorted(original_adata.obs[groupby2].unique())
        default_color_map2 = dict(zip(all_unique_labels_secondary, tab20_colors[:len(all_unique_labels_secondary)]))
        
        if groupby2_label_color_map:
            color_map2 = groupby2_label_color_map
        else:
            color_map2 = default_color_map2

    # Setup subplot configuration
    y_position = list(range(len(valid_genes)))
    heatmap_height = 40 * len(valid_genes)
    bar_chart_height1 = 30
    bar_chart_height2 = 30 if has_secondary else 0
    
    if has_secondary:
        total_height = [heatmap_height, bar_chart_height1, bar_chart_height2]
        rows = 3
    else:
        total_height = [heatmap_height, bar_chart_height1]
        rows = 2
    
    total_x_range = sum(value_list1)

    # Create subplot
    fig = make_subplots(
        rows=rows, cols=1,
        row_heights=total_height,
        shared_xaxes=True,
        vertical_spacing=0.02
    )

    # Set colorbar range based on transformation type
    if z_score:
        z_max = max(abs(heatmap_gene_matrix.min()), abs(heatmap_gene_matrix.max()))
        zmin, zmax = -z_max, z_max
        zmid = 0
    else:
        zmin = heatmap_gene_matrix.min()
        zmax = heatmap_gene_matrix.max() 
        zmid = None

    # Create colorbar configuration
    colorbar_config = dict(
        title=f'Expression(log)' if log else f'Expression(z-score)' if z_score else 'Expression',
        len=0.4 if has_secondary else 0.5,
        y=1,
        yanchor='top'
    )

    # Add main heatmap
    heatmap = go.Heatmap(
        z=heatmap_gene_matrix,
        x=list(range(len(sorted_heatmap_df))),
        y=valid_genes,
        colorscale=color_map,
        colorbar=colorbar_config,
        hoverinfo='skip',
        zmin=zmin,
        zmax=zmax,
        zmid=zmid,
    )
    fig.add_trace(heatmap, row=1, col=1)

    # Add boundary lines if requested
    if boundary is not False:
        boundary_width = boundary if isinstance(boundary, (int, float)) else 1
        boundary_positions = np.cumsum(value_list1[:-1]).tolist()
        for pos in boundary_positions:
            fig.add_shape(
                type="line",
                x0=pos, y0=-0.5, x1=pos, y1=len(valid_genes) - 0.5,
                xref='x1', yref='y1',
                line=dict(color="rgba(0, 0, 0, 0.8)", width=boundary_width),
                row=1, col=1
            )

    # Add primary annotation bar
    for i, label in enumerate(label_list1):
        fig.add_trace(go.Bar(
            x=[value_list1[i]],
            marker_color=color_map1[label],
            name=f"{label}",
            hovertemplate=f'<b>{label}</b><br>Count: {value_list1[i]}<extra></extra>',
            orientation='h',
            showlegend=False,
            width=1,
        ), row=2, col=1)

    # Add secondary annotation bar if present
    if has_secondary:
        secondary_labels_shown = set()
        for i, label in enumerate(label_list2):
            secondary_label = label2_dict[label]
            show_legend = secondary_label not in secondary_labels_shown
            if show_legend:
                secondary_labels_shown.add(secondary_label)

            fig.add_trace(go.Bar(
                x=[value_list2[i]],
                marker_color=color_map2[secondary_label],
                name=f"{secondary_label}",
                hovertemplate=f'<b>{secondary_label}</b><br>Count: {value_list2[i]}<extra></extra>',
                orientation='h',
                showlegend=False,
                width=2,
            ), row=3, col=1)

    # Create legend annotations
    legend_annotations = []
    legend_start_x = 1.01
    legend_start_y = 0.5
    current_y = legend_start_y

    # Primary annotation legend
    legend_annotations.append(dict(
        x=legend_start_x, y=current_y,
        xref="paper", yref="paper",
        text=f"<b>{groupby1}</b>",
        showarrow=False,
        font=dict(size=12, color='black'),
        xanchor='left', yanchor='top'
    ))
    
    current_y -= 0.08
    for label in label_list1:
        legend_annotations.append(dict(
            x=legend_start_x, y=current_y,
            xref="paper", yref="paper",
            text=f"<span style='color:{color_map1[label]}'>■</span> {label}",
            showarrow=False,
            font=dict(size=12),
            xanchor='left', yanchor='middle'
        ))
        current_y -= 0.04

    # Secondary annotation legend
    if has_secondary:
        current_y -= 0.04  # Add extra spacing between sections
        legend_annotations.append(dict(
            x=legend_start_x, y=current_y,
            xref="paper", yref="paper",
            text=f"<b>{groupby2}</b>",
            showarrow=False,
            font=dict(size=12, color='black'),
            xanchor='left', yanchor='top'
        ))
        
        current_y -= 0.08
        unique_secondary_labels_sorted = sorted(set(label2_dict.values()))
        for label in unique_secondary_labels_sorted:
            legend_annotations.append(dict(
                x=legend_start_x, y=current_y,
                xref="paper", yref="paper",
                text=f"<span style='color:{color_map2[label]}'>■</span> {label}",
                showarrow=False,
                font=dict(size=12),
                xanchor='left', yanchor='middle'
            ))
            current_y -= 0.04

    # Update layout
    layout_updates = {
        'barmode': 'stack',
        'showlegend': False,
        'plot_bgcolor': 'white',
        'paper_bgcolor': 'white',
        'annotations': legend_annotations,
        'hovermode': 'closest',
        'height': max(450, sum(total_height)),
        'margin': dict(t=50, b=150 if not has_secondary else 150, l=50, r=200),
        'xaxis': dict(
            range=[0, total_x_range],
            constrain='domain',
            constraintoward='center'
        ),
        'yaxis': dict(
            tickmode='array',
            tickvals=y_position,
            ticktext=valid_genes,
            showgrid=False,
            zeroline=False,
            tickfont=dict(size=14),
            constrain='domain',
            constraintoward='middle'
        ),
        'yaxis2': dict(
            tickmode='array',
            tickvals=[0],
            ticktext=[f"<b>{groupby1}</b>"],
            showgrid=False,
            zeroline=False,
            tickfont=dict(size=14),
            constrain='domain',
            constraintoward='middle'
        )
    }
    
    # Add secondary y-axis configuration if needed
    if has_secondary:
        layout_updates['xaxis2'] = dict(
            range=[0, total_x_range],  
            visible=False,
            constrain='domain',
            constraintoward='center'
        )
        layout_updates['xaxis3'] = dict(
            range=[0, total_x_range],  
            visible=False,
            constrain='domain',
            constraintoward='center'
        )
        layout_updates['yaxis3'] = dict(
            tickmode='array',
            tickvals=[0],
            ticktext=[f"<b>{groupby2}</b>"],
            showgrid=False,
            zeroline=False,
            tickfont=dict(size=14),
            constrain='domain',
            constraintoward='middle'
        )
    else:
        layout_updates['xaxis2'] = dict(
            range=[0, total_x_range], 
            visible=False,
            constrain='domain',
            constraintoward='center'
        )

    fig.update_layout(**layout_updates)
    fig.update_xaxes(tickangle=90, visible=False)

    return fig

def plot_heatmap2_continuous(
    adata, genes, groupby1, continuous_key, labels=None, log=False, z_score=False, 
    color_map='Viridis', groupby1_label_color_map=None, max_cells=50000, n_bins=10000
):
    is_backed = hasattr(adata, "isbacked") and adata.isbacked
    filtered_obs = adata.obs
    filtered_obs_names = adata.obs_names

    if labels:
        cell_indices = adata.obs[groupby1].isin(labels)
        if is_backed:
            cell_indices_array = np.where(cell_indices)[0]
            filtered_obs = adata.obs.iloc[cell_indices_array]
            filtered_obs_names = adata.obs_names[cell_indices_array]
        else:
            adata = adata[cell_indices].copy()
            filtered_obs = adata.obs
            filtered_obs_names = adata.obs_names

    valid_genes = [gene for gene in genes if gene in adata.var_names]
    if not valid_genes:
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

    if labels and is_backed:
        gene_df_list = []
        for gene in valid_genes:
            gene_expr = extract_gene_expression(adata, gene)
            gene_expr_filtered = gene_expr[cell_indices_array]
            gene_df_list.append(pd.Series(gene_expr_filtered, name=gene, index=filtered_obs_names))
        gene_df = pd.concat(gene_df_list, axis=1)
    else:
        gene_df = extract_multiple_genes(adata, valid_genes)

    # 转换
    if log:
        for gene in valid_genes:
            gene_df[gene] = apply_transformation(gene_df[gene], method='log1p')
    if z_score:
        for gene in valid_genes:
            gene_df[gene] = apply_transformation(gene_df[gene], method='z_score')

    gene_df.insert(0, 'CellID', filtered_obs_names)
    label_df = pd.DataFrame(filtered_obs[[groupby1, continuous_key]])
    label_df.insert(0, 'CellID', filtered_obs_names)
    heatmap_df = pd.merge(gene_df, label_df, on='CellID')

    sorted_heatmap_df = heatmap_df.sort_values(continuous_key)
    use_binning = len(sorted_heatmap_df) > max_cells
    if use_binning:
        sorted_heatmap_df = bin_cells_for_heatmap(
            sorted_heatmap_df, valid_genes, groupby1, n_bins, continuous_key=continuous_key
        )
    heatmap_gene_matrix = sorted_heatmap_df[valid_genes].values.T

    # Get ALL unique labels from original dataset for consistent color assignment
    # Use the original adata (before any filtering) to get all labels
    all_unique_labels = sorted(adata.obs[groupby1].unique())
    
    # Import the consistent color config
    from guanaco.data_loader import color_config
    
    # Create color map using ALL labels (not just filtered ones)
    default_color_map1 = {
        label: color_config[i % len(color_config)] for i, label in enumerate(all_unique_labels)
    }
    
    # Use provided color map or default
    if groupby1_label_color_map:
        color_map1 = groupby1_label_color_map
    else:
        color_map1 = default_color_map1

    heatmap_height = 40 * len(valid_genes)
    bar_chart_height = 30
    continuous_bar_height = 30

    fig = make_subplots(
        rows=3, cols=1,
        row_heights=[heatmap_height, continuous_bar_height, bar_chart_height],
        shared_xaxes=True,
        vertical_spacing=0.01
    )

    # Set colorbar range based on transformation type
    if z_score:
        # For z-score, use symmetric range centered at 0
        z_max = max(abs(heatmap_gene_matrix.min()), abs(heatmap_gene_matrix.max()))
        zmin, zmax = -z_max, z_max
        zmid = 0
    else:
        # For other transformations, use data range
        zmin = heatmap_gene_matrix.min()
        zmax = heatmap_gene_matrix.max() 
        zmid = None

    heatmap = go.Heatmap(
        z=heatmap_gene_matrix,
        x=list(range(len(sorted_heatmap_df))),
        y=valid_genes,
        colorscale=color_map,
        colorbar=dict(
            title=f'Expression(log)' if log else f'Expression(z-score)' if z_score else 'Expression',
            len=0.4,
            y=1,
            yanchor='top'
        ),
        hoverinfo='skip',
        zmin=zmin,
        zmax=zmax,
        zmid=zmid,
    )
    fig.add_trace(heatmap, row=1, col=1)

    continuous_values = sorted_heatmap_df[continuous_key].values.reshape(1, -1)
    fig.add_trace(go.Heatmap(
        z=continuous_values,
        x=list(range(len(sorted_heatmap_df))),
        y=[continuous_key],
        colorscale='Viridis',
        showscale=False,
        hovertemplate=f'{continuous_key}: %{{z:.4f}}<extra></extra>'
    ), row=2, col=1)

    category_to_num = {cat: i for i, cat in enumerate(sorted_heatmap_df[groupby1].unique())}
    category_values = sorted_heatmap_df[groupby1].map(category_to_num).values.reshape(1, -1)

    unique_categories = sorted(sorted_heatmap_df[groupby1].unique())
    n_categories = len(unique_categories)
    colorscale = []
    for i, cat in enumerate(unique_categories):
        colorscale.append([i/n_categories, color_map1[cat]])
        colorscale.append([(i+1)/n_categories, color_map1[cat]])

    fig.add_trace(go.Heatmap(
        z=category_values,
        x=list(range(len(sorted_heatmap_df))),
        y=[groupby1],
        colorscale=colorscale,
        showscale=False,
        text=[[cat for cat in sorted_heatmap_df[groupby1]]],
        hovertemplate='%{text}<extra></extra>',
        zmin=-0.5,
        zmax=n_categories-0.5
    ), row=3, col=1)

    legend_annotations = []
    legend_start_x = 1.01
    legend_start_y = 0.5
    legend_annotations.append(dict(
        x=legend_start_x, y=legend_start_y,
        xref="paper", yref="paper",
        text=f"<b>{groupby1}</b>",
        showarrow=False,
        font=dict(size=12, color='black'),
        xanchor='left', yanchor='top'
    ))
    for i, label in enumerate(unique_labels):
        y_pos = legend_start_y - 0.08 - (i * 0.04)
        legend_annotations.append(dict(
            x=legend_start_x, y=y_pos,
            xref="paper", yref="paper",
            text=f"<span style='color:{color_map1[label]}'>■</span> {label}",
            showarrow=False,
            font=dict(size=14),
            xanchor='left', yanchor='middle'
        ))

    fig.update_layout(
        barmode='stack',
        showlegend=False,
        plot_bgcolor='white',
        paper_bgcolor='white',
        xaxis=dict(visible=False, constrain='domain', constraintoward='center'),
        xaxis2=dict(visible=False, constrain='domain', constraintoward='center'),
        xaxis3=dict(visible=False, constrain='domain', constraintoward='center'),
        yaxis=dict(
            tickmode='array',
            tickvals=list(range(len(valid_genes))),
            ticktext=valid_genes,
            showgrid=False,
            zeroline=False,
            tickfont=dict(size=14),
            constrain='domain',
            constraintoward='middle'
        ),
        yaxis2=dict(
            tickmode='array',
            tickvals=[0],
            ticktext=[continuous_key],
            showgrid=False,
            zeroline=False,
            tickfont=dict(size=12),
            constrain='domain',
            constraintoward='middle'
        ),
        yaxis3=dict(
            tickmode='array',
            tickvals=[0],
            ticktext=[groupby1],
            showgrid=False,
            zeroline=False,
            tickfont=dict(size=12),
            constrain='domain',
            constraintoward='middle'
        ),
        annotations=legend_annotations,
        hovermode='closest', 
        height=max(450, sum([heatmap_height, continuous_bar_height, bar_chart_height])),
        margin=dict(t=50, b=50, l=50, r=150),
    )
    return fig

