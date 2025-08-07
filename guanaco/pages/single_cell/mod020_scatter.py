"""
Scatter plot layout module for single cell visualization
"""
import json
from pathlib import Path
from dash import dcc, html
import dash_bootstrap_components as dbc
import plotly.express as px

# Import configs
from guanaco.pages.single_cell.config import scatter_config

# Load color palettes
cvd_color_path = Path(__file__).parent / "cvd_color.json"
with open(cvd_color_path, "r") as f:
    palette_json = json.load(f)
palette_names = list(palette_json["color_palettes"].keys())


def initialize_scatter_components(adata):
    """Initialize scatter plot components based on available embeddings"""
    embedding_prefixes = {
        "X_umap": "UMAP", "X_pca": "PCA", "X_tsne": "t-SNE", "X_diffmap": "DiffMap",
        "X_phate": "PHATE", "X_draw_graph_fa": "FA"
    }
    
    obsm_list = list(adata.obsm.keys())
    embedding_columns = {
        key: [f'{embedding_prefixes.get(key, key.upper())}{i+1}' for i in range(adata.obsm[key].shape[1])]
        for key in obsm_list
    }
    default_embedding = obsm_list[-1]
    default_columns = embedding_columns[default_embedding]
    
    return obsm_list, embedding_columns, default_embedding, default_columns


def create_control_components(adata, prefix):
    """Create control components for dimension reduction selection"""
    obsm_list, _, default_embedding, default_columns = initialize_scatter_components(adata)
    clustering_dropdown = dcc.Dropdown(
        id=f'{prefix}-clustering-dropdown',
        options=[{'label': key.replace("X_", "").upper(), 'value': key} for key in obsm_list],
        value='X_umap' if 'X_umap' in obsm_list else default_embedding,
        placeholder="Select Clustering Method",
        style={'marginBottom': '15px','fontSize': '14px'},
        clearable=False
    )
    
    coordinates_dropdowns = html.Div([
        html.Label("X-axis:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
        dcc.Dropdown(id=f'{prefix}-x-axis', options=[{'label': col, 'value': col} for col in default_columns], value=default_columns[0]),
        html.Label("Y-axis:", style={'fontWeight': 'bold', 'marginTop': '15px', 'marginBottom': '5px'}),
        dcc.Dropdown(id=f'{prefix}-y-axis', options=[{'label': col, 'value': col} for col in default_columns], value=default_columns[1])
    ], style={'marginBottom': '15px','fontSize': '14px'})
    
    return clustering_dropdown, coordinates_dropdowns


def generate_annotation_dropdown(anno_list, prefix):
    """Generate annotation dropdown for scatter plot"""
    return dcc.Dropdown(id=f'{prefix}-annotation-dropdown', 
    options=[{'label': label, 'value': label} for label in anno_list],
    placeholder="Search annotations or genes...", 
    value = anno_list[0] if anno_list else None,
    style={'marginBottom': '15px'})


def generate_scatter_gene_selection(init_gene_list, prefix):
    """Generate gene selection dropdown for scatter plot"""
    # Find the first gene in the list (skipping annotations which come first)
    default_value = None
    for item in init_gene_list:
        # Assuming genes start after the first 20 or so items (which are annotations)
        # We can check if item looks like a gene name (often uppercase or mixed case)
        # For safety, just use item at index 20 if list is long enough
        if len(init_gene_list) > 20:
            default_value = init_gene_list[20]  # First gene after annotations
            break
    if default_value is None:
        default_value = init_gene_list[0] if init_gene_list else None
    
    return dcc.Dropdown(id=f'{prefix}-scatter-gene-selection', options=[{'label': label, 'value': label} for label in init_gene_list], 
    value = default_value, placeholder="Search and select a gene...", style={'marginBottom': '15px'})


def create_global_metadata_filter(adata, prefix):
    """Create global metadata filter component"""
    
    # Get all categorical metadata columns
    categorical_columns = []
    for col in adata.obs.columns:
        if adata.obs[col].dtype == 'category' or adata.obs[col].dtype == 'object':
            unique_vals = adata.obs[col].unique()
            if len(unique_vals) < 100:  # Reasonable limit for categorical
                categorical_columns.append(col)
    
    # Create filter components
    filter_components = []
    
    for col in categorical_columns:
        unique_values = sorted([str(val) for val in adata.obs[col].unique()])
        
        filter_component = html.Div([
            html.Label(f"{col}:", style={'fontWeight': 'bold', 'marginRight': '10px', 'minWidth': '120px'}),
            dcc.Dropdown(
                id={'type': f'{prefix}-global-metadata-filter', 'column': col},
                options=[{'label': val, 'value': val} for val in unique_values],
                value=unique_values,  # All selected by default
                multi=True,
                placeholder=f"Select {col}...",
                style={'flex': '1', 'minWidth': '200px'}
            )
        ], style={
            'display': 'flex',
            'alignItems': 'center',
            'marginBottom': '8px',
            'padding': '5px'
        })
        
        filter_components.append(filter_component)
    
    # Global filter panel
    global_filter_panel = html.Div([
        # Header with cell count
        html.Div([
            html.H5("🔍 Global Data Filter", 
                   style={'margin': '0', 'color': '#2c3e50', 'display': 'inline-block'}),
            html.Div([
                html.Span("Active cells: ", style={'fontWeight': 'bold'}),
                html.Span(
                    id=f'{prefix}-global-cell-count',
                    children=f"{adata.n_obs:,}",
                    style={'color': '#27ae60', 'fontWeight': 'bold', 'fontSize': '16px'}
                ),
                html.Span(f" / {adata.n_obs:,} total", style={'color': '#7f8c8d'}),
                html.Span(
                    id=f'{prefix}-filter-preview',
                    children="",
                    style={'marginLeft': '15px', 'color': '#e67e22', 'fontSize': '14px', 'fontStyle': 'italic'}
                )
            ], style={'display': 'inline-block', 'marginLeft': '20px'})
        ], style={'marginBottom': '15px'}),
        
        # Collapsible filter section
        dbc.Collapse([
            html.Div(filter_components, style={'maxHeight': '300px', 'overflowY': 'auto'}),
            
            # Quick action buttons
            html.Div([
                dbc.Button("Select All", id=f'{prefix}-select-all-filters', 
                          color="success", size="sm", style={'marginRight': '10px'}),
                dbc.Button("Clear All", id=f'{prefix}-clear-all-filters', 
                          color="warning", size="sm", style={'marginRight': '10px'}),
                dbc.Button("Apply Filter", id=f'{prefix}-apply-global-filter', 
                          color="primary", size="sm", 
                          style={'fontWeight': 'bold', 'minWidth': '100px'})
            ], style={'textAlign': 'center', 'marginTop': '15px'})
            
        ], id=f'{prefix}-global-filter-collapse', is_open=False),
        
        # Toggle button
        html.Div([
            dbc.Button(
                "▼ Show Filters",
                id=f'{prefix}-toggle-global-filter',
                color="link",
                size="sm",
                style={'padding': '5px 10px', 'textDecoration': 'none'}
            )
        ], style={'textAlign': 'center', 'marginTop': '10px'}),
        
        # Hidden store for filtered data (start empty for better performance)
        dcc.Store(id=f'{prefix}-global-filtered-data', data={'cell_indices': None, 'n_cells': adata.n_obs})
        
    ], style={
        'backgroundColor': '#f8f9fa',
        'border': '2px solid #dee2e6',
        'borderRadius': '10px',
        'padding': '15px',
        'marginBottom': '20px',
        'boxShadow': '0 2px 4px rgba(0,0,0,0.1)'
    })
    
    return global_filter_panel


def generate_scatter_layout(adata, prefix):
    """Generate the complete scatter plot layout"""
    
    scatter_transformation_selection = html.Div([
        dbc.RadioItems(
            id=f'{prefix}-scatter-log-or-zscore',
            options=[
                {'label': 'None', 'value':None},
                {'label': 'Log', 'value': 'log'}
            ],
            value=None,
            inline=True,
            style={'fontSize': '14px'} 
        )
    ])

    scatter_order_selection = html.Div([
        dbc.RadioItems(
            id=f'{prefix}-plot-order',
            options=[
                {'label': 'Max-1st', 'value': 'max'},
                {'label': 'Min-1st', 'value': 'min'},
                {'label': 'Original', 'value': 'original'},
                {'label': 'Random', 'value': 'random'},
            ],
            value='max',
            inline=True,
            style={'fontSize': '14px'}
        )
    ])

    colorscales = px.colors.named_colorscales()
    # Dropdown for color map selection
    color_map_dropdown = dcc.Dropdown(
        id=f"{prefix}-scatter-color-map-dropdown",
        options=[{"label": scale, "value": scale} for scale in colorscales],
        value="viridis",  # Default colorscale
        style={"marginBottom": "10px"},
        clearable = False
    )

    color_map_discrete_dropdown = dcc.Dropdown(
        id=f"{prefix}-discrete-color-map-dropdown",
        options=[{"label": name, "value": name} for name in palette_names],
        value=None,
        placeholder="Default color",
        style={"marginBottom": "10px"},
        clearable=True,
        className="custom-dropdown",
    )

    # Slider for marker size
    marker_size_slider = dcc.Slider(
        id=f'{prefix}-marker-size-slider',
        min=1,
        max=10,
        value=5,  # Default marker size
        marks = None,
        tooltip={"placement": "bottom", "always_visible": True},
        className="dbc-slider"
    )

    # Slider for opacity
    opacity_slider = dcc.Slider(
        id=f"{prefix}-opacity-slider",
        min=0.1,
        max=1.0,
        value=1,  # Default opacity
        marks = None,
        tooltip={"placement": "bottom", "always_visible": True},
        className="dbc-slider"
    )
    
    # toggle for scatter plot legend
    scatter_legend_toggle = dbc.RadioItems(
        id=f'{prefix}-scatter-legend-toggle',
        options=[
            {'label': 'Show', 'value': 'show'},
            {'label': 'Hide', 'value': 'hide'}
        ],
        value='hide',
        inline=True,
        style={'fontSize': '14px'} 
    )
    
    # Toggle for axis show
    axis_toggle = dbc.RadioItems(
        id=f'{prefix}-axis-toggle',
        options=[
            {'label': 'Show', 'value': True},
            {'label': 'Hide', 'value': False}
        ],
        value=True,
        inline=True,
        style={'fontSize': '14px'}
    )

    graphic_control = html.Div(
        id=f"{prefix}-controls-container",
        children=[
            html.Div(
                [
                    html.Label("Plot order: ", className="control-label"),
                    scatter_order_selection,
                ],
                style={'marginBottom': '15px'}
            ),
            html.Div(
                [
                    html.Label("Continous ColorMap:",  className="control-label"),
                    color_map_dropdown,
                ],
                className="dbc",
                style={'marginBottom': '15px'}
            ),
            html.Div(
                [
                    html.Label("Discrete ColorMap:",  className="control-label"),
                    color_map_discrete_dropdown,
                ],
                className="dbc",
                style={'marginBottom': '15px'}
            ),
            html.Div(
                [
                    html.Label("Marker Size:",  className="control-label"),
                    marker_size_slider,
                ],
                style={'marginBottom': '15px'}
            ),
            html.Div(
                [
                    html.Label("Opacity:",  className="control-label"),
                    opacity_slider,
                ],
                style={'marginBottom': '15px'}
            ),
            html.Div(
                [
                    html.Label("Legend on data:",  className="control-label"),
                    scatter_legend_toggle,
                ],
                style={'marginBottom': '15px'}
            ),
            html.Div(
                [
                    html.Label("Axis:",  className="control-label"),
                    axis_toggle,
                ],
                style={'marginBottom': '15px'}
            ),
        ],
        style={'display': 'none'}  # hidden by default
    )
    
    anno_list = adata.obs.columns.tolist()
    sample_genes = adata.var_names[:20].tolist()
    anno_list.extend(sample_genes)

    clustering_dropdown, coordinates_dropdowns = create_control_components(adata, prefix)
    
    layout = html.Div([
        # Add Store component to hold selected cells data
        dcc.Store(id=f'{prefix}-selected-cells-store'),
        
        # Add global metadata filter at the top
        create_global_metadata_filter(adata, prefix),
        
        dbc.Row([
            # Left column: Controls
            dbc.Col(
            html.Div(
                [
                    html.Div(
                        [
                            html.Label("Dimension Reduction:", className = "control-label"),
                            clustering_dropdown,
                        ],
                        className="dbc",
                        style={'marginBottom': '15px'}
                    ),
                    html.Div(
                        [
                            html.Div(coordinates_dropdowns, id=f'{prefix}-coordinates-dropdowns')
                        ],
                        className="dbc",
                        style={'marginBottom': '15px'}
                    ),
                    html.Div(
                        [
                            html.Label("Transformation:",  className = "control-label"),
                            scatter_transformation_selection,
                        ],
                        style={'marginBottom': '15px'}
                    ),
                    html.Button("More controls", id=f"{prefix}-toggle-button", n_clicks=0, style={'marginBottom': '10px','border':'1px solid','borderRadius': '5px'}),

                    graphic_control,
                ],
            ),
            xs=12, sm=12, md=3, lg=3, xl=2,  # Reduced width for controls
            style={"borderRight": "1px solid #ddd", "padding": "10px"}
        ),

        # Combined scatter plot column
        dbc.Col(
            html.Div(
                [
                    # Controls row for both plots
                    dbc.Row([
                        dbc.Col([
                            html.Label("Search Annotation/Gene:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                            generate_annotation_dropdown(anno_list=anno_list, prefix=prefix),
                        ], width=6),
                        dbc.Col([
                            html.Label("Search Annotation/Gene:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                            generate_scatter_gene_selection(init_gene_list=anno_list, prefix=prefix),
                        ], width=6),
                    ], style={'marginBottom': '10px'}),
                    
                    # Co-expression controls
                    dbc.Row([
                        dbc.Col([
                            # Placeholder to maintain alignment
                        ], width=6),
                        dbc.Col([
                            # Add toggle for co-expression mode
                            dbc.RadioItems(
                                id=f'{prefix}-coexpression-toggle',
                                options=[
                                    {'label': 'Single Gene', 'value': 'single'},
                                    {'label': 'Co-expression', 'value': 'coexpression'}
                                ],
                                value='single',
                                inline=True,
                                style={'marginBottom': '10px', 'fontSize': '14px'}
                            ),
                            # Second gene selection (hidden by default)
                            html.Div(
                                [
                                    html.Label("Second Gene:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                                    dcc.Dropdown(
                                        id=f'{prefix}-scatter-gene2-selection',
                                        options=[{'label': label, 'value': label} for label in anno_list],
                                        value=anno_list[21] if len(anno_list) > 21 else (anno_list[1] if len(anno_list) > 1 else None),
                                        placeholder="Search and select second gene...",
                                        style={'marginBottom': '10px'}
                                    ),
                                ],
                                id=f'{prefix}-gene2-container',
                                style={'display': 'none'}
                            ),
                            # Threshold sliders 
                            html.Div(
                                [
                                    html.Label("Expression Thresholds:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                                    html.Div([
                                        html.Label("Gene 1 Threshold:", style={'fontSize': '12px'}),
                                        dcc.Slider(
                                            id=f'{prefix}-gene1-threshold-slider',
                                            min=0,
                                            max=1,
                                            value=0.5,
                                            marks=None,
                                            tooltip={"placement": "bottom", "always_visible": True},
                                            className="dbc-slider"
                                        ),
                                    ], style={'marginBottom': '10px'}),
                                    html.Div([
                                        html.Label("Gene 2 Threshold:", style={'fontSize': '12px'}),
                                        dcc.Slider(
                                            id=f'{prefix}-gene2-threshold-slider',
                                            min=0,
                                            max=1,
                                            value=0.5,
                                            marks=None,
                                            tooltip={"placement": "bottom", "always_visible": True},
                                            className="dbc-slider"
                                        ),
                                    ], style={'marginBottom': '10px'}),
                                ],
                                id=f'{prefix}-threshold-container',
                                style={'display': 'none'}
                            ),
                        ], width=6),
                    ]),
                    
                    # Combined scatter plot and buttons container using flexbox
                    html.Div([
                        # Plot container that grows to fill available space
                        html.Div([
                            dcc.Loading(
                                id=f"{prefix}-loading-combined-scatter",
                                type="circle",
                                children=dcc.Graph(
                                    id=f'{prefix}-combined-scatter', 
                                    config=scatter_config,
                                    style={
                                        "height": "100%",  # Fill container height
                                        "width": "100%",   # Full width of container
                                    },
                                    responsive=True  # Make plot responsive to container size
                                ),
                                style={"height": "100%", "width": "100%"},
                            ),
                        ], style={
                            "flex": "1",  # Grow to fill available space
                            "minHeight": "500px",  # Minimum height
                            "height": "calc(100vh - 350px)",  # Dynamic height based on viewport
                            "overflow": "hidden",  # Prevent overflow
                            "position": "relative"  # For proper positioning
                        }),
                        
                        # Keep the original annotation and gene scatter plots as hidden for backward compatibility
                        html.Div([
                            dcc.Graph(id=f'{prefix}-annotation-scatter', style={'display': 'none'}),
                            dcc.Graph(id=f'{prefix}-gene-scatter', style={'display': 'none'}),
                        ]),
                        
                        # Buttons directly below the plot (no extra spacing)
                        html.Div([
                            dbc.Row([
                                dbc.Col([
                                    dbc.Button(
                                        "Update other Plots",
                                        id=f"{prefix}-update-plots-button",
                                        color="primary",
                                        n_clicks=0,
                                        style={'width': '100%'}
                                    ),
                                ], xs=12, sm=12, md=8, lg=8, xl=8),
                                dbc.Col([
                                    dbc.InputGroup([
                                        dbc.DropdownMenu(
                                            [
                                                dbc.DropdownMenuItem("Cell IDs (.txt)", id=f"{prefix}-download-cellids"),
                                                dbc.DropdownMenuItem("Subset AnnData (.h5ad)", id=f"{prefix}-download-adata"),
                                            ],
                                            label="Download",
                                            color="secondary",
                                            id=f"{prefix}-download-menu",
                                            disabled=True
                                        ),
                                        dcc.Download(id=f"{prefix}-download-cells-data")
                                    ], style={'width': '100%'})
                                ], xs=12, sm=12, md=4, lg=4, xl=4)
                            ], style={'marginTop': '10px'}),
                            html.Div(id=f"{prefix}-selection-status", style={'textAlign': 'center', 'marginTop': '5px'})
                        ], style={
                            "flex": "0 0 auto",  # Don't grow, stay at natural size
                            "marginTop": "-10px"  # Negative margin to pull buttons closer to legend
                        }),
                    ], style={
                        "display": "flex",
                        "flexDirection": "column",
                        "height": "auto",  # Let content determine height
                        "width": "100%",  # Use full width available
                    }),
                ],
                className="dbc",
                style={'marginBottom': '20px'}
            ),
            xs=12, sm=12, md=9, lg=9, xl=10  # Increased width for scatter plot
        ),
    ])
    ])

    return layout