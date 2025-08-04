import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
import plotly.express as px
from guanaco.pages.single_cell.cellplotly.gene_extraction_utils import (
    extract_multiple_genes, apply_transformation, bin_cells_for_heatmap
)




def plot_heatmap2_continuous(adata, genes, groupby1, continuous_key, labels=None, log=False, z_score=False, 
                           color_map='Viridis', groupby1_label_color_map=None, max_cells=50000, n_bins=10000):
    """
    Simplified heatmap for continuous secondary annotation (e.g., pseudotime).
    Orders cells by the continuous value and shows primary annotation as colored bars below.
    """
    # Filter data based on selected labels
    if labels:
        cell_indices = adata.obs[groupby1].isin(labels)
        adata = adata[cell_indices]
    
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
 
    # Use centralized gene extraction
    gene_df = extract_multiple_genes(adata, valid_genes)
    
    # Apply transformations using centralized utility
    if log:
        for gene in valid_genes:
            gene_df[gene] = apply_transformation(gene_df[gene], method='log1p')
    if z_score:
        for gene in valid_genes:
            gene_df[gene] = apply_transformation(gene_df[gene], method='z_score')
    
    # Add cell IDs
    gene_df.insert(0, 'CellID', adata.obs_names)

    # Process labels
    label_df = pd.DataFrame(adata.obs[[groupby1, continuous_key]])
    label_df.insert(0, 'CellID', adata.obs_names)
    
    heatmap_df = pd.merge(gene_df, label_df, on='CellID')
    
    # Sort by continuous value (e.g., pseudotime) and apply binning if needed
    sorted_heatmap_df = heatmap_df.sort_values(continuous_key)
    
    # Apply binning if dataset is large
    use_binning = len(sorted_heatmap_df) > max_cells
    if use_binning:
        sorted_heatmap_df = bin_cells_for_heatmap(sorted_heatmap_df, valid_genes, groupby1, n_bins, continuous_key=continuous_key)
    
    heatmap_gene_matrix = sorted_heatmap_df[valid_genes].values.T
    
    # Get unique labels for coloring
    unique_labels = sorted_heatmap_df[groupby1].unique()
    colors1 = px.colors.qualitative.Plotly
    color_map1 = dict(zip(unique_labels, colors1))
    
    if groupby1_label_color_map:
        color_map1 = groupby1_label_color_map
    
    # Calculate heights
    heatmap_height = 40 * len(valid_genes)
    bar_chart_height = 30
    continuous_bar_height = 30
    
    # Create subplots (3 rows: heatmap, continuous values, categorical annotation)
    fig = make_subplots(
        rows=3, cols=1,
        row_heights=[heatmap_height, continuous_bar_height, bar_chart_height],
        shared_xaxes=True,
        vertical_spacing=0.01
    )
    
    
    # Add heatmap
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
        hoverinfo='skip',  # Disable hover
        zmin=heatmap_gene_matrix.min(),
        zmax=heatmap_gene_matrix.max(),
    )
    fig.add_trace(heatmap, row=1, col=1)
    
    # Add continuous values as a heatmap strip
    continuous_values = sorted_heatmap_df[continuous_key].values.reshape(1, -1)
    fig.add_trace(go.Heatmap(
        z=continuous_values,
        x=list(range(len(sorted_heatmap_df))),
        y=[continuous_key],
        colorscale='Viridis',
        showscale=False,
        hovertemplate=f'{continuous_key}: %{{z:.4f}}<extra></extra>'
    ), row=2, col=1)
    
    # Add categorical annotation as a colored heatmap strip
    # Map categories to numeric values for visualization
    category_to_num = {cat: i for i, cat in enumerate(sorted_heatmap_df[groupby1].unique())}
    category_values = sorted_heatmap_df[groupby1].map(category_to_num).values.reshape(1, -1)
    
    # Create custom colorscale based on the color mapping
    unique_categories = sorted(sorted_heatmap_df[groupby1].unique())
    n_categories = len(unique_categories)
    
    # Create discrete colorscale
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
    
    # Create legend annotations
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
    
    # Add legend items
    for i, label in enumerate(unique_labels):
        y_pos = legend_start_y - 0.08 - (i * 0.04)
        legend_annotations.append(dict(
            x=legend_start_x, y=y_pos,
            xref="paper", yref="paper",
            text=f"<span style='color:{color_map1[label]}'>■</span> {label}",
            showarrow=False,
            font=dict(size=10),
            xanchor='left', yanchor='middle'
        ))
    
    # Update layout
    fig.update_layout(
        barmode='stack',
        showlegend=False,
        plot_bgcolor='white',
        paper_bgcolor='white',
        xaxis=dict(
            visible=False,
            constrain='domain',
            constraintoward='center'
        ),
        xaxis2=dict(
            visible=False,
            constrain='domain',
            constraintoward='center'
        ),
        xaxis3=dict(
            visible=False,
            constrain='domain',
            constraintoward='center'
        ),
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


def plot_heatmap2(adata, genes, groupby1, groupby2, labels=None, log=False, z_score=False,boundary=False,color_map='Viridis',groupby1_label_color_map=None,groupby2_label_color_map=None):
    
    # Filter data based on selected labels
    if labels:
        cell_indices = adata.obs[groupby1].isin(labels)
        adata = adata[cell_indices]
    
    # Helper function to check if annotation is continuous
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
    
    # If secondary annotation is continuous, use the specialized function
    if is_continuous_annotation(adata, groupby2):
        return plot_heatmap2_continuous(adata, genes, groupby1, groupby2, labels, log, z_score, 
                                      color_map, groupby1_label_color_map)
    
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
            gene_expression_matrix = adata.X[:, gene_indices].toarray()
        else:
            gene_expression_matrix = adata.X[:, gene_indices]
    else:
        # Original code for non-backed AnnData
        adata_selected = adata[:, valid_genes]
        gene_expression_matrix = adata_selected.X.toarray()
    
    if log:
        gene_expression_matrix = np.log1p(gene_expression_matrix)  # log1p for log(1 + x) transformation

    if z_score:
        gene_expression_matrix = (gene_expression_matrix - gene_expression_matrix.mean(axis=0)) / gene_expression_matrix.std(axis=0)

    gene_expression_df = pd.DataFrame(gene_expression_matrix, columns=valid_genes)
    gene_expression_df.insert(0, 'CellID', adata.obs_names)

    # Process labels and combine groupby1 and groupby2
    label_df = pd.DataFrame(adata.obs[[groupby1, groupby2]])
    label_df.insert(0, 'CellID', adata.obs_names)
    
    heatmap_df = pd.merge(gene_expression_df, label_df, on='CellID')
    
    # Sort by groupby1 using the input labels order, then by groupby2
    if labels:
        # Create categorical with specified order
        heatmap_df[groupby1] = pd.Categorical(heatmap_df[groupby1], categories=labels, ordered=True)
    
    # Sort by both annotations to ensure consistent ordering
    sorted_heatmap_df = heatmap_df.sort_values([groupby1, groupby2])
    
    # Combine groupby1 and groupby2 labels
    sorted_heatmap_df['combined'] = sorted_heatmap_df[groupby1].astype(str) + '_' + sorted_heatmap_df[groupby2].astype(str)
    
    heatmap_gene_matrix = sorted_heatmap_df[valid_genes].values.T
    
    # Calculate groupby labels distribution
    label_list1 = sorted_heatmap_df[groupby1].unique().tolist()
    
    value_list1 = []
    for item in label_list1:
        value_list1.append(sorted_heatmap_df[sorted_heatmap_df[groupby1] == item].shape[0])
    
    label_list2 = sorted_heatmap_df['combined'].unique().tolist()
    
    label2_dict = {item: item.split('_')[-1] for item in label_list2}
    
    value_list2 = []
    for item in label_list2:
        value_list2.append(sorted_heatmap_df[sorted_heatmap_df['combined'] == item].shape[0])

    # Use different Plotly color palettes for distinct themes
    colors1 = px.colors.qualitative.Plotly 
    colors2 = px.colors.qualitative.Set1 

    color_map1 = dict(zip(label_list1, colors1))
    color_map2 = dict(zip(label2_dict.values(), colors2))

    if groupby1_label_color_map:
        color_map1 = groupby1_label_color_map
    if groupby2_label_color_map:
        color_map2 = groupby2_label_color_map
        
    # y_position and y_labels setup
    y_position = list(range(len(valid_genes)))
    y_labels = valid_genes
    
    # Calculate the height dynamically - remove annotation rows
    heatmap_height = 40 * len(valid_genes)
    bar_chart_height1 = 30
    bar_chart_height2 = 30
    
    # Create a list for row_heights directly (only 3 rows now)
    total_height = [heatmap_height, bar_chart_height1, bar_chart_height2]

    # Create subplots
    fig = make_subplots(
        rows=3, cols=1,
        row_heights=total_height,
        shared_xaxes=True,
        vertical_spacing=0.02
    )

    total_x_range = sum(value_list1)


    # Add heatmap
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
        hoverinfo='skip',  # Disable hover
        zmin=heatmap_gene_matrix.min(),
        zmax=heatmap_gene_matrix.max(),
    )
    fig.add_trace(heatmap, row=1, col=1)

    if boundary is not False:
        boundary_width = boundary
        boundary_positions = np.cumsum(value_list1[:-1]).tolist()
        for pos in boundary_positions:
            fig.add_shape(
                type="line",
                x0=pos, y0=-0.5, x1=pos, y1=len(valid_genes) - 0.5,
                xref='x1', yref='y1',
                line=dict(color="rgba(0, 0, 0, 0.8)", width=boundary_width),
                row=1, col=1
            )

    x_pos = 0
    for i, label in enumerate(label_list1):
        fig.add_trace(go.Bar(
            x=[value_list1[i]],
            marker_color=color_map1[label],
            name=f"{label}",
            hovertext=f'{value_list1[i]} ({label})',
            hoverinfo='text',
            orientation='h',
            showlegend=False,
            width=1,
        ), row=2, col=1)
        x_pos += value_list1[i]

    x2_pos = 0
    # Create a set to track unique secondary labels for legend
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
            hovertext=f'{value_list2[i]} ({secondary_label})',
            hoverinfo='text',
            orientation='h',
            showlegend=False,
            width=2,
        ), row=3, col=1)

        # Remove individual labels under bars - legend will show them
        x2_pos += value_list2[i]

            
    # Create custom legend using annotations
    legend_annotations = []
    
    # Calculate positions for two-column legend matching heatmap height
    legend_start_x = 1.01
    legend_start_y = 0.5  # Top of heatmap area
    
    # Column positions with more space between them
    col1_x = legend_start_x
    col2_x = legend_start_x + 0.14  # Adjusted column spacing
    
    # Add first group header
    legend_annotations.append(dict(
        x=col1_x, y=legend_start_y,
        xref="paper", yref="paper",
        text=f"<b>{groupby1}</b>",
        showarrow=False,
        font=dict(size=12, color='black'),
        xanchor='left', yanchor='top'
    ))
    
    # Add first group items with proper spacing below header
    for i, label in enumerate(label_list1):
        y_pos = legend_start_y - 0.08 - (i * 0.04)  # Start lower, smaller spacing
        legend_annotations.append(dict(
            x=col1_x, y=y_pos,
            xref="paper", yref="paper",
            text=f"<span style='color:{color_map1[label]}'>■</span> {label}",
            showarrow=False,
            font=dict(size=10),
            xanchor='left', yanchor='middle'
        ))
    
    # Add second group header
    legend_annotations.append(dict(
        x=col2_x, y=legend_start_y,
        xref="paper", yref="paper",
        text=f"<b>{groupby2}</b>",
        showarrow=False,
        font=dict(size=12, color='black'),
        xanchor='left', yanchor='top'
    ))
    
    # Add second group items with proper spacing below header
    unique_secondary_labels = sorted(list(secondary_labels_shown))
    for i, label in enumerate(unique_secondary_labels):
        y_pos = legend_start_y - 0.08 - (i * 0.04)  # Start lower, smaller spacing
        legend_annotations.append(dict(
            x=col2_x, y=y_pos,
            xref="paper", yref="paper",
            text=f"<span style='color:{color_map2[label]}'>■</span> {label}",
            showarrow=False,
            font=dict(size=10),
            xanchor='left', yanchor='middle'
        ))

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
        xaxis3=dict(
            range=[0, total_x_range],  
            visible=False,
            constrain='domain',
            constraintoward='center'
        ),
        yaxis=dict(
            tickmode='array',
            tickvals=y_position,
            ticktext=y_labels,
            showgrid=False,  
            zeroline=False,  
            tickfont=dict(size=14),  
            titlefont=dict(size=14),
            constrain='domain',
            constraintoward='middle'
        ),
        yaxis2=dict(
            tickmode='array',
            tickvals=[0],
            ticktext=[f"<b>{groupby1}</b>"],  
            showgrid=False,
            zeroline=False,
            tickfont=dict(size=14),
            constrain='domain',
            constraintoward='middle'
        ),
        yaxis3=dict(
            tickmode='array',
            tickvals=[0],
            ticktext=[f"<b>{groupby2}</b>"],  
            showgrid=False,
            zeroline=False,
            tickfont=dict(size=14),
            constrain='domain',
            constraintoward='middle'
        ),
        annotations=legend_annotations,
        hovermode='closest',  
        height=max(450, sum(total_height)),
        margin=dict(t=50, b=150, l=50, r=200),  # Further increased right margin for legend
    )
    fig.update_xaxes(tickangle=90, visible=False)

    return fig