import plotly.graph_objs as go
import pandas as pd
import numpy as np
from plotly.subplots import make_subplots
from .gene_extraction_utils import extract_gene_expression, apply_transformation
from guanaco.data_loader import color_config
import json
from pathlib import Path

def plot_combined_scatter_subplots(
    adata, adata_full, embedding_key, left_feature, right_feature, 
    x_axis=None, y_axis=None,
    transformation=None, order=None, 
    left_type='auto', right_type='auto', right_plot_mode='single',
    second_gene=None, threshold1=0.5, threshold2=0.5,
    color_map='Viridis', color_map_discrete=None,
    marker_size=5, opacity=1,
    legend_show='hide', axis_show=True
):

    # Helper function to determine feature type
    def determine_type(feature, type_hint='auto'):
        if type_hint != 'auto':
            return type_hint
        if feature in adata.obs.columns:
            # Check if it's categorical
            if (adata.obs[feature].dtype == 'category' or 
                adata.obs[feature].dtype == 'object' or
                adata.obs[feature].nunique() < 50):
                return 'categorical'
            else:
                return 'continuous'
        elif feature in adata.var_names:
            return 'gene'
        else:
            return 'continuous'  # Default
    
    # Determine actual types for both features
    left_type_actual = determine_type(left_feature, left_type)
    right_type_actual = determine_type(right_feature, right_type)
    

    # Create subplots with shared axes and fixed column widths
    subplot_title_right = f"<b>{right_feature}</b>"
    if right_type_actual == 'gene' and right_plot_mode == 'coexpression' and second_gene:
        subplot_title_right = f"<b>Co-expression: {right_feature} & {second_gene}</b>"
    
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=[
            f"<b>{left_feature}</b>",
            subplot_title_right
        ],
        shared_xaxes=True,
        shared_yaxes=True,
        horizontal_spacing=0.15,
        column_widths=[0.5, 0.5]  # Equal width for both plots
    )
    
    # Manually position subplot titles to be very close to plots
    for i, title_annotation in enumerate(fig.layout.annotations):

        # Position titles over respective subplots
        if i == 0:  # Left subplot title
            x_pos = 0.225  # Center of left subplot (0 to 0.45, center = 0.225)
        else:  # Right subplot title
            x_pos = 0.775  # Center of right subplot (0.55 to 1, center = 0.775)
        
        title_annotation.update(
            font=dict(size=18),
            xanchor="center", 
            yanchor="bottom",
            x=x_pos, 
            y=0.98
        )
    
    # Get embedding data
    embedding_prefixes = {
        'X_umap': 'UMAP', 'X_pca': 'PCA', 'X_tsne': 't-SNE',
        'X_diffmap': 'DiffMap', 'X_phate': 'PHATE', 'X_draw_graph_fa': 'FA'
    }
    embedding_prefix = embedding_prefixes.get(embedding_key, embedding_key.upper())
    embedding_data = adata.obsm[embedding_key]
    
    # Set column names for the embedding
    num_dimensions = embedding_data.shape[1]
    embedding_columns = [f'{embedding_prefix}{i + 1}' for i in range(num_dimensions)]
    embedding_df = pd.DataFrame(embedding_data, columns=embedding_columns)
    
    # Default x and y axis
    x_axis = x_axis or embedding_columns[0]
    y_axis = y_axis or (embedding_columns[1] if len(embedding_columns) > 1 else embedding_columns[0])
    
    # Helper function to add a plot to the figure
    def add_plot(feature, feature_type, row, col):
        if feature_type == 'categorical':
            # Categorical plot
            on_data = legend_show == 'show'  # Show legend on data when 'show'
            embedding_df[feature] = adata.obs[feature].values
            
            # Get all unique labels from FULL dataset for consistent color mapping
            if feature in adata_full.obs.columns:
                all_unique_labels = sorted(adata_full.obs[feature].unique())
            else:
                all_unique_labels = sorted(embedding_df[feature].unique())
            
            # Get labels present in current (filtered) data
            unique_labels_filtered = sorted(embedding_df[feature].unique())
            
            
            if color_map_discrete:
                cvd_color_path = Path(__file__).parent.parent / "cvd_color.json"
                with open(cvd_color_path, "r") as f:
                    palette_json = json.load(f)
                color_palette = palette_json["color_palettes"][color_map_discrete]
            else:
                # Use color_config for consistency with violin and heatmap plots
                color_palette = color_config
            
            # Create consistent color mapping based on ALL labels from full dataset
            label_to_color = {
                label: color_palette[i % len(color_palette)]
                for i, label in enumerate(all_unique_labels)
            }
            
            # Add grey background (non-selectable)
            fig.add_trace(go.Scattergl(
                x=embedding_df[x_axis],
                y=embedding_df[y_axis],
                mode='markers',
                marker=dict(
                    size=marker_size,
                    color='lightgrey',
                    opacity=opacity * 0.3,
                ),
                name='Background',
                hoverinfo='skip',
                showlegend=False,
                visible=True,
                selectedpoints=[],  # Disable selection for background
                unselected=dict(marker=dict(opacity=0))  # Make unselected invisible
            ), row=row, col=col)
            
            # Add traces for each category present in filtered data
            for label in unique_labels_filtered:
                mask = embedding_df[feature] == label
                masked_indices = np.arange(len(adata))[mask]
                
                # Legend logic:
                # - Always show the right-side legend (regardless of on_data setting)
                # - For left plot (col==1), always show legend in the right side
                # - For right plot (col==2), show legend only if it's a different categorical feature
                if col == 1:
                    show_in_legend = True  # Always show legend for left plot
                else:
                    # Show legend for right plot if it's categorical and different from left
                    show_in_legend = (feature != left_feature)
                
                # Use flat array for customdata (Plotly may not serialize 2D arrays properly)
                customdata_array = masked_indices.tolist()
                
                fig.add_trace(go.Scattergl(
                    x=embedding_df.loc[mask, x_axis],
                    y=embedding_df.loc[mask, y_axis],
                    mode='markers',
                    marker=dict(
                        size=marker_size,
                        color=label_to_color[label],
                        opacity=opacity,
                    ),
                    name=str(label),
                    customdata=customdata_array,
                    hoverinfo='skip',
                    showlegend=show_in_legend,
                    legendgroup=str(label) if col == 1 else str(label) + '_right',
                    selectedpoints=None,
                    selected=dict(marker=dict(opacity=1)),
                    unselected=dict(marker=dict(opacity=0.2))
                ), row=row, col=col)
            
            # Add text annotations if legend_show == 'on data'
            if on_data:
                for label in unique_labels_filtered:
                    mask = embedding_df[feature] == label
                    # Calculate centroid of the cluster
                    centroid_x = embedding_df.loc[mask, x_axis].median()
                    centroid_y = embedding_df.loc[mask, y_axis].median()
                    
                    # Add annotation at centroid
                    fig.add_annotation(
                        x=centroid_x,
                        y=centroid_y,
                        text=f"<b>{label}</b>",
                        showarrow=False,
                        font=dict(size=12, color='black'),
                        xanchor='center',
                        yanchor='middle',
                        row=row,
                        col=col
                    )
        else:
            # Continuous or gene plot
            if feature in adata.obs.columns and feature_type == 'continuous':
                values = adata.obs[feature].values
            else:
                values = extract_gene_expression(adata, feature)
            
            if transformation:
                values = apply_transformation(values, transformation, copy=False)
            
            embedding_df[feature] = values
            embedding_df['_cell_idx'] = np.arange(len(embedding_df))
            
            # Sort order
            if order == 'max':
                embedding_df_sorted = embedding_df.sort_values(by=feature)
            elif order == 'min':
                embedding_df_sorted = embedding_df.sort_values(by=feature, ascending=False)
            elif order == 'random':
                embedding_df_sorted = embedding_df.sample(frac=1, random_state=315).reset_index(drop=True)
            else:
                embedding_df_sorted = embedding_df
            
            # Configure colorbar position based on plot location

            if col == 1:
                colorbar_config = dict(
                    title=f"{feature}<br>{transformation if transformation else ''}",
                    orientation='h',
                    len=0.35,
                    thickness=15,
                    x=0.225,
                    xanchor='center',
                    y=-0.18,
                    yanchor='top'
                )
            else:
                colorbar_config = dict(
                    title=f"{feature}<br>{transformation if transformation else ''}",
                    orientation='v',
                    len=0.5,
                    thickness=15,
                    x=1.02,
                    xanchor='left',
                    y=0.5,
                    yanchor='middle'
                )
            
            fig.add_trace(go.Scattergl(
                x=embedding_df_sorted[x_axis],
                y=embedding_df_sorted[y_axis],
                mode='markers',
                marker=dict(
                    color=embedding_df_sorted[feature],
                    colorscale=color_map,
                    cmin=embedding_df_sorted[feature].min(),
                    cmax=embedding_df_sorted[feature].max(),
                    size=marker_size,
                    opacity=opacity,
                    colorbar=colorbar_config
                ),
                customdata=embedding_df_sorted['_cell_idx'].values.tolist(),
                hoverinfo='skip',
                showlegend=False,
                selectedpoints=None,
                selected=dict(marker=dict(opacity=1)),
                unselected=dict(marker=dict(opacity=0.2))
            ), row=row, col=col)
    

    # Dual plot mode
    # LEFT PLOT
    add_plot(left_feature, left_type_actual, 1, 1)
    
    # RIGHT PLOT - handle special case of coexpression
    if right_type_actual == 'gene' and right_plot_mode == 'coexpression' and second_gene:
        # Co-expression plot
        gene1_expr = extract_gene_expression(adata, right_feature)
        gene2_expr = extract_gene_expression(adata, second_gene)
        
        if transformation:
            gene1_expr = apply_transformation(gene1_expr, transformation, copy=True)
            gene2_expr = apply_transformation(gene2_expr, transformation, copy=True)
        
        # Categorize cells
        gene1_expressed = gene1_expr > threshold1
        gene2_expressed = gene2_expr > threshold2
        
        categories = np.zeros(len(adata), dtype='object')
        categories[:] = 'Neither'
        categories[gene1_expressed & ~gene2_expressed] = f'{right_feature} only'
        categories[~gene1_expressed & gene2_expressed] = f'{second_gene} only'
        categories[gene1_expressed & gene2_expressed] = 'Co-expressed'
        
        embedding_df['category'] = categories
        
        # Color mapping
        coexpr_colors = {
            'Neither': '#E8E8E8',
            f'{right_feature} only': '#648fff',
            f'{second_gene} only': '#ffb000',
            'Co-expressed': '#dc267f'
        }
        
        category_order = ['Neither', f'{right_feature} only', f'{second_gene} only', 'Co-expressed']
        
        for category in category_order:
            mask = embedding_df['category'] == category
            if mask.any():
                fig.add_trace(go.Scattergl(
                    x=embedding_df.loc[mask, x_axis],
                    y=embedding_df.loc[mask, y_axis],
                    mode='markers',
                    marker=dict(
                        size=marker_size,
                        color=coexpr_colors.get(category, '#808080'),
                        opacity=opacity,
                    ),
                    name=category,
                    customdata=np.arange(len(adata))[mask].tolist(),
                    hoverinfo='skip',
                    showlegend=False,
                    selectedpoints=None,
                    selected=dict(marker=dict(opacity=1)),
                    unselected=dict(marker=dict(opacity=0.2))
                ), row=1, col=2)
    else:
        # Regular plot for right side
        add_plot(right_feature, right_type_actual, 1, 2)
    
    # Configure legend
    # Check if we need legends (either plot is categorical)
    # We always have 2 columns in the subplot layout
    has_categorical = left_type_actual == 'categorical' or right_type_actual == 'categorical'
    
    # Always show the right-side legend when there are categorical features
    if has_categorical:
        # If both are categorical and different, we need to handle dual legends
        if (left_type_actual == 'categorical' and 
            right_type_actual == 'categorical' and left_feature != right_feature):
            # Position legend to accommodate both
            legend_config = dict(
                orientation='h',
                itemsizing='constant',
                x=0, y=-0.19,  # Positioned below x-axis without overlapping
                xanchor='left',
                yanchor='top',
                bgcolor='rgba(0,0,0,0)',
                bordercolor='rgba(0,0,0,0)',
                borderwidth=0,
                itemclick='toggle',
                itemdoubleclick='toggleothers',
                font=dict(size=10),
                xref='paper',
                yref='paper',
                traceorder='grouped',  # Group legends by trace group
                tracegroupgap=20,  # Add gap between groups
                itemwidth=30
            )
        else:
            # Single legend configuration
            legend_config = dict(
                orientation='h',
                itemsizing='constant',
                x=0, y=-0.19,  # Positioned below x-axis without overlapping
                xanchor='left',
                yanchor='top',
                bgcolor='rgba(0,0,0,0)',
                bordercolor='rgba(0,0,0,0)',
                borderwidth=0,
                itemclick='toggle',
                itemdoubleclick='toggleothers',
                font=dict(size=10),
                xref='paper',
                yref='paper',
                traceorder='normal',
                tracegroupgap=5,
                itemwidth=30
            )
    else:
        legend_config = None
    

    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        height=None,  # Let container control height
        autosize=True,
        showlegend=has_categorical,  # Always show legend when categorical features present
        legend=legend_config,
        uirevision='constant',
        margin=dict(t=50, b=50, l=60, r=60),  # Reduced bottom margin
        dragmode='zoom'
    )
    
    # Calculate data ranges for synchronization
    x_min = embedding_df[x_axis].min()
    x_max = embedding_df[x_axis].max()
    y_min = embedding_df[y_axis].min()
    y_max = embedding_df[y_axis].max()
    
    # Add padding to ranges
    x_padding = (x_max - x_min) * 0.05
    y_padding = (y_max - y_min) * 0.05
    
    x_range = [x_min - x_padding, x_max + x_padding]
    y_range = [y_min - y_padding, y_max + y_padding]
    
    # Update axes for both subplots with synchronized ranges
    fig.update_xaxes(
        title=x_axis,
        showgrid=False,
        zeroline=False,
        range=x_range,
        matches='x2',
        constrain='domain',
        showline=True,
        linewidth=2,
        linecolor='black',
        tickfont=dict(color='black' if axis_show else 'rgba(0,0,0,0)'),
        row=1, col=1
    )
    
    fig.update_xaxes(
        title='',
        showgrid=False,
        zeroline=False,
        range=x_range,
        matches='x',
        constrain='domain',
        showline=True,
        linewidth=2,
        linecolor='black',
        tickfont=dict(color='black' if axis_show else 'rgba(0,0,0,0)'),
        row=1, col=2
    )
    
    fig.update_yaxes(
        title=y_axis,
        showgrid=False,
        zeroline=False,
        range=y_range,
        matches='y2',
        scaleanchor='x',
        scaleratio=1,
        constrain='domain',
        showline=True,
        linewidth=2,
        linecolor='black',
        tickfont=dict(color='black' if axis_show else 'rgba(0,0,0,0)'),
        row=1, col=1
    )
    
    fig.update_yaxes(
        title='',
        showgrid=False,
        zeroline=False,
        range=y_range,
        matches='y',
        scaleanchor='x2',
        scaleratio=1,
        constrain='domain',
        showline=True,
        linewidth=2,
        linecolor='black',
        tickfont=dict(color='black' if axis_show else 'rgba(0,0,0,0)'),
        row=1, col=2
    )
    
    # Override domain for better use of space
    fig.update_layout(
        xaxis=dict(domain=[0, 0.45]),      # Left plot: 0% to 45%
        xaxis2=dict(domain=[0.55, 1.0]),   # Right plot: 55% to 100% (10% gap)
    )
    
    return fig