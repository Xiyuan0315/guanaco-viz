import plotly.graph_objs as go
import pandas as pd
import numpy as np
import plotly.express as px
from plotly.subplots import make_subplots
from .gene_extraction_utils import extract_gene_expression, apply_transformation


def plot_combined_scatter_subplots(
    adata, adata_full, embedding_key, annotation, gene, 
    x_axis=None, y_axis=None,
    transformation=None, order=None, 
    annotation_type='auto', gene_plot_type='single',
    gene2=None, threshold1=0.5, threshold2=0.5,
    color_map='Viridis', color_map_discrete=None,
    marker_size=5, opacity=1, 
    legend_show='right', axis_show=True
):
    """
    Create combined subplot with annotation and gene scatter plots that share axes.
    
    Parameters:
    -----------
    adata : AnnData
        The annotated data object
    adata_full : AnnData
        The full dataset for consistent color mapping
    embedding_key : str
        Key for embedding in adata.obsm
    annotation : str
        Column name from adata.obs or gene name for left plot
    gene : str
        Gene name for right plot
    annotation_type : str
        'continuous', 'categorical', or 'auto' (auto-detect)
    gene_plot_type : str
        'single' for single gene or 'coexpression' for two genes
    gene2 : str
        Second gene name (only used if gene_plot_type='coexpression')
    Other parameters match the individual plotting functions
    """
    
    # Determine plot types
    if annotation_type == 'auto':
        if annotation in adata.obs.columns:
            # Check if it's categorical
            if (adata.obs[annotation].dtype == 'category' or 
                adata.obs[annotation].dtype == 'object' or
                adata.obs[annotation].nunique() < 50):
                annotation_type = 'categorical'
            else:
                annotation_type = 'continuous'
        elif annotation in adata.var_names:
            annotation_type = 'gene'
        else:
            annotation_type = 'continuous'  # Default
    
    # Check if we're creating a single plot (same annotation and gene)
    single_plot_mode = (annotation == gene)
    
    if single_plot_mode:
        # Create a single centered plot instead of subplots
        fig = make_subplots(
            rows=1, cols=1,
            subplot_titles=[f"<b>{annotation}</b>"]
        )
        # Use row=1, col=1 for all traces in single plot mode
        target_row, target_col = 1, 1
    else:
        # Create subplots with proper spacing
        fig = make_subplots(
            rows=1, cols=2,
            shared_xaxes='all',  # Use 'all' for stronger sharing
            shared_yaxes='all',  # Use 'all' for stronger sharing
            subplot_titles=[f"<b>{annotation}</b>", f"<b>{gene}</b>" if gene_plot_type == 'single' else f"<b>Co-expression: {gene} & {gene2}</b>"],
            horizontal_spacing=0.20,  # Increase spacing to 20% to prevent overlap
            column_widths=[0.40, 0.40]  # Reduce plot widths to 40% each
        )
    
    # Position subplot titles properly
    for i, title_annotation in enumerate(fig.layout.annotations):
        if single_plot_mode:
            # Center title for single plot
            title_annotation.update(
                font=dict(size=18),
                xanchor="center", 
                yanchor="bottom",
                x=0.5, 
                y=0.98
            )
        else:
            # Position titles for dual mode with new domains
            if i == 0:  # Left subplot title
                x_pos = 0.20  # Center of left plot (0 to 0.40)
            else:  # Right subplot title  
                x_pos = 0.80  # Center of right plot (0.60 to 1.0)
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
    
    # Track all traces for selection
    all_customdata = []
    
    if single_plot_mode:
        # SINGLE PLOT MODE: Use only the annotation/gene for plotting
        if annotation_type == 'categorical':
            # Categorical single plot
            on_data = legend_show == 'on data'
            embedding_df[annotation] = adata.obs[annotation].values
            
            # Get all unique labels from FULL dataset for consistent color mapping
            if annotation in adata_full.obs.columns:
                all_unique_labels = sorted(adata_full.obs[annotation].unique())
            else:
                all_unique_labels = sorted(embedding_df[annotation].unique())
            
            # Get labels present in current (filtered) data
            unique_labels_filtered = sorted(embedding_df[annotation].unique())
            
            # Load color palettes for consistent color assignment
            import json
            from pathlib import Path
            cvd_color_path = Path(__file__).parent.parent / "cvd_color.json"
            with open(cvd_color_path, "r") as f:
                palette_json = json.load(f)
            
            if color_map_discrete:
                color_palette = palette_json["color_palettes"][color_map_discrete]
            else:
                color_palette = palette_json["color_palettes"]["Okabe_Ito"]
            
            # Create consistent color mapping based on ALL labels from full dataset
            label_to_color = {
                label: color_palette[i % len(color_palette)]
                for i, label in enumerate(all_unique_labels)
            }
            
            # Add grey background
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
                visible=True
            ), row=target_row, col=target_col)
            
            # Add traces for each category present in filtered data
            for label in unique_labels_filtered:
                mask = embedding_df[annotation] == label
                cell_indices = np.arange(len(adata))[mask]
                all_customdata.extend(cell_indices)
                
                show_in_legend = not on_data
                
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
                    customdata=cell_indices,
                    hoverinfo='skip',
                    showlegend=show_in_legend,
                    legendgroup=str(label),
                    selectedpoints=None,
                    selected=dict(marker=dict(opacity=1)),
                    unselected=dict(marker=dict(opacity=0.2))
                ), row=target_row, col=target_col)
        
        else:
            # Continuous single plot (gene or continuous annotation)
            if annotation in adata.obs.columns:
                values = adata.obs[annotation].values
            else:
                values = extract_gene_expression(adata, annotation)
            
            if transformation:
                values = apply_transformation(values, transformation, copy=False)
            
            embedding_df[annotation] = values
            embedding_df['_cell_idx'] = np.arange(len(embedding_df))
            all_customdata = embedding_df['_cell_idx'].values
            
            # Sort order
            if order == 'max':
                embedding_df_sorted = embedding_df.sort_values(by=annotation)
            elif order == 'min':
                embedding_df_sorted = embedding_df.sort_values(by=annotation, ascending=False)
            elif order == 'random':
                embedding_df_sorted = embedding_df.sample(frac=1, random_state=315).reset_index(drop=True)
            else:
                embedding_df_sorted = embedding_df
            
            fig.add_trace(go.Scattergl(
                x=embedding_df_sorted[x_axis],
                y=embedding_df_sorted[y_axis],
                mode='markers',
                marker=dict(
                    color=embedding_df_sorted[annotation],
                    colorscale=color_map,
                    cmin=embedding_df_sorted[annotation].min(),
                    cmax=embedding_df_sorted[annotation].max(),
                    size=marker_size,
                    opacity=opacity,
                    colorbar=dict(
                        title=f"{annotation}<br>{transformation if transformation else ''}",
                        orientation='h',  # Horizontal colorbar
                        len=0.6,
                        thickness=15,
                        x=0.5,  # Center of plot
                        xanchor='center',
                        y=-0.10,
                        yanchor='top'
                    )
                ),
                customdata=embedding_df_sorted['_cell_idx'],
                hoverinfo='skip',
                selectedpoints=None,
                selected=dict(marker=dict(opacity=1)),
                unselected=dict(marker=dict(opacity=0.2))
            ), row=target_row, col=target_col)
    else:
        # DUAL PLOT MODE: Original subplot logic (annotation on left, gene on right)
        # LEFT PLOT: Annotation (same as before)
        if annotation_type == 'categorical':
            # Use existing categorical plotting logic
            on_data = legend_show == 'on data'
            embedding_df[annotation] = adata.obs[annotation].values
            
            # Get all unique labels from FULL dataset for consistent color mapping
            if annotation in adata_full.obs.columns:
                all_unique_labels = sorted(adata_full.obs[annotation].unique())
            else:
                all_unique_labels = sorted(embedding_df[annotation].unique())
            
            # Get labels present in current (filtered) data
            unique_labels_filtered = sorted(embedding_df[annotation].unique())
            
            # Load color palettes for consistent color assignment
            import json
            from pathlib import Path
            cvd_color_path = Path(__file__).parent.parent / "cvd_color.json"
            with open(cvd_color_path, "r") as f:
                palette_json = json.load(f)
            
            if color_map_discrete:
                color_palette = palette_json["color_palettes"][color_map_discrete]
            else:
                color_palette = palette_json["color_palettes"]["Okabe_Ito"]
            
            # Create consistent color mapping based on ALL labels from full dataset
            label_to_color = {
                label: color_palette[i % len(color_palette)]
                for i, label in enumerate(all_unique_labels)
            }
            
            # Add grey background
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
            ), row=1, col=1)
            
            # Add traces for each category present in filtered data
            for label in unique_labels_filtered:
                mask = embedding_df[annotation] == label
                cell_indices = np.arange(len(adata))[mask]
                all_customdata.extend(cell_indices)
                
                show_in_legend = not on_data
                
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
                    customdata=cell_indices,
                    hoverinfo='skip',
                    showlegend=show_in_legend,
                    legendgroup=str(label),
                    selectedpoints=None,
                    selected=dict(marker=dict(opacity=1)),
                    unselected=dict(marker=dict(opacity=0.2)),
                    xaxis='x',
                    yaxis='y'
                ), row=1, col=1)
        
        else:  # continuous or gene
            # Extract values
            if annotation in adata.obs.columns:
                values = adata.obs[annotation].values
            else:
                values = extract_gene_expression(adata, annotation)
            
            if transformation:
                values = apply_transformation(values, transformation, copy=False)
            
            embedding_df[annotation] = values
            embedding_df['_cell_idx'] = np.arange(len(embedding_df))
            all_customdata = embedding_df['_cell_idx'].values
            
            # Sort order
            if order == 'max':
                embedding_df_sorted = embedding_df.sort_values(by=annotation)
            elif order == 'min':
                embedding_df_sorted = embedding_df.sort_values(by=annotation, ascending=False)
            elif order == 'random':
                embedding_df_sorted = embedding_df.sample(frac=1, random_state=315).reset_index(drop=True)
            else:
                embedding_df_sorted = embedding_df
            
            fig.add_trace(go.Scattergl(
                x=embedding_df_sorted[x_axis],
                y=embedding_df_sorted[y_axis],
                mode='markers',
                marker=dict(
                    color=embedding_df_sorted[annotation],
                    colorscale=color_map,
                    cmin=embedding_df_sorted[annotation].min(),
                    cmax=embedding_df_sorted[annotation].max(),
                    size=marker_size,
                    opacity=opacity,
                    colorbar=dict(
                        title=f"{annotation}<br>{transformation if transformation else ''}",
                        orientation='v',
                        len=0.5,
                        thickness=15,
                        x=-0.08,
                        xanchor='right',
                        y=0.5,
                        yanchor='middle'
                    )
                ),
                customdata=embedding_df_sorted['_cell_idx'],
                hoverinfo='skip',
                selectedpoints=None,
                selected=dict(marker=dict(opacity=1)),
                unselected=dict(marker=dict(opacity=0.2)),
            ), row=1, col=1)
        
        # RIGHT PLOT: Gene (same as before but always goes to col=2)
        if gene_plot_type == 'coexpression' and gene2:
            # Co-expression plot
            gene1_expr = extract_gene_expression(adata, gene)
            gene2_expr = extract_gene_expression(adata, gene2)
            
            if transformation:
                gene1_expr = apply_transformation(gene1_expr, transformation, copy=True)
                gene2_expr = apply_transformation(gene2_expr, transformation, copy=True)
            
            # Categorize cells
            gene1_expressed = gene1_expr > threshold1
            gene2_expressed = gene2_expr > threshold2
            
            categories = np.zeros(len(adata), dtype='object')
            categories[:] = 'Neither'
            categories[gene1_expressed & ~gene2_expressed] = f'{gene} only'
            categories[~gene1_expressed & gene2_expressed] = f'{gene2} only'
            categories[gene1_expressed & gene2_expressed] = 'Co-expressed'
            
            embedding_df['category'] = categories
            
            # Color mapping
            coexpr_colors = {
                'Neither': '#E8E8E8',
                f'{gene} only': '#648fff',
                f'{gene2} only': '#ffb000',
                'Co-expressed': '#dc267f'
            }
            
            category_order = ['Neither', f'{gene} only', f'{gene2} only', 'Co-expressed']
            
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
                        customdata=np.arange(len(adata))[mask],
                        hoverinfo='skip',
                        selectedpoints=None,
                        selected=dict(marker=dict(opacity=1)),
                        unselected=dict(marker=dict(opacity=0.2))
                    ), row=1, col=2)
        
        else:
            # Single gene plot
            gene_expr = extract_gene_expression(adata, gene)
            
            if transformation:
                gene_expr = apply_transformation(gene_expr, transformation, copy=False)
            
            embedding_df[gene] = gene_expr
            embedding_df['_cell_idx'] = np.arange(len(embedding_df))
            
            # Sort order
            if order == 'max':
                embedding_df_sorted = embedding_df.sort_values(by=gene)
            elif order == 'min':
                embedding_df_sorted = embedding_df.sort_values(by=gene, ascending=False)
            elif order == 'random':
                embedding_df_sorted = embedding_df.sample(frac=1, random_state=315).reset_index(drop=True)
            else:
                embedding_df_sorted = embedding_df
            
            fig.add_trace(go.Scattergl(
                x=embedding_df_sorted[x_axis],
                y=embedding_df_sorted[y_axis],
                mode='markers',
                marker=dict(
                    color=embedding_df_sorted[gene],
                    colorscale=color_map,
                    cmin=embedding_df_sorted[gene].min(),
                    cmax=embedding_df_sorted[gene].max(),
                    size=marker_size,
                    opacity=opacity,
                    colorbar=dict(
                        title=f"{gene}<br>{transformation if transformation else ''}",
                        orientation='v',  # Vertical colorbar
                        len=0.8,
                        thickness=15,
                        x=1.02,
                        xanchor='left',
                        y=0.5,
                        yanchor='middle'
                    )
                ),
                customdata=embedding_df_sorted['_cell_idx'],
                hoverinfo='skip',
                showlegend=False,
                selectedpoints=None,
                selected=dict(marker=dict(opacity=1)),
                unselected=dict(marker=dict(opacity=0.2))
            ), row=1, col=2)
    
    # Configure legend
    if legend_show == 'right' and annotation_type == 'categorical':
        legend_config = dict(
            orientation='h',
            itemsizing='constant',
            x=0, y=-0.25,  # Increased margin to prevent button overlap
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
    
    # Update layout
    if single_plot_mode:
        fig.update_layout(
            plot_bgcolor='white',
            paper_bgcolor='white',
            height=None,  # Let container control height
            autosize=True,
            showlegend=(legend_show == 'right' and annotation_type == 'categorical'),
            legend=legend_config,
            margin=dict(t=20, b=120),  # Increased bottom margin for legend spacing
            uirevision='constant',
            dragmode='zoom',  # Enable zoom mode to match scatter_config
            xaxis=dict(
                title=x_axis,
                showgrid=False,
                zeroline=False,
                constrain='domain',
                showline=True,
                linewidth=2,
                linecolor='black',
                tickfont=dict(color='black' if axis_show else 'rgba(0,0,0,0)'),
                scaleanchor='y',
                scaleratio=1
            ),
            yaxis=dict(
                title=y_axis,
                showgrid=False,
                zeroline=False,
                constrain='domain',
                showline=True,
                linewidth=2,
                linecolor='black',
                tickfont=dict(color='black' if axis_show else 'rgba(0,0,0,0)')
            )
        )
    else:
        # Subplot layout configuration
        fig.update_layout(
            plot_bgcolor='white',
            paper_bgcolor='white',
            height=None,  # Let container control height
            autosize=True,
            showlegend=(legend_show == 'right' and annotation_type == 'categorical'),
            legend=legend_config,
            uirevision='constant',
            margin=dict(t=50, b=100, l=50, r=50),
            dragmode='zoom',
            # Set domains to prevent overlap while maintaining axis sharing
            xaxis=dict(
                domain=[0, 0.40],
                matches='x2',  # Explicitly match x2 axis
                showgrid=False,
                zeroline=False,
                constrain='domain',
                showline=True,
                linewidth=2,
                linecolor='black',
                tickfont=dict(color='black' if axis_show else 'rgba(0,0,0,0)')
            ),
            xaxis2=dict(
                domain=[0.60, 1.0],
                showgrid=False,
                zeroline=False,
                constrain='domain',
                showline=True,
                linewidth=2,
                linecolor='black',
                tickfont=dict(color='black' if axis_show else 'rgba(0,0,0,0)')
            ),
            yaxis=dict(
                domain=[0, 1],
                matches='y2',  # Explicitly match y2 axis
                showgrid=False,
                zeroline=False,
                constrain='domain',
                showline=True,
                linewidth=2,
                linecolor='black',
                tickfont=dict(color='black' if axis_show else 'rgba(0,0,0,0)')
            ),
            yaxis2=dict(
                domain=[0, 1],
                showgrid=False,
                zeroline=False,
                constrain='domain',
                showline=True,
                linewidth=2,
                linecolor='black',
                tickfont=dict(color='black' if axis_show else 'rgba(0,0,0,0)')
            )
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
        
        # Update titles and ranges for both axes
        fig.update_xaxes(title=x_axis, range=x_range)
        fig.update_yaxes(title=y_axis, range=y_range)
    
    # Axis synchronization is handled by scaleanchor properties in the layout
    
    return fig
