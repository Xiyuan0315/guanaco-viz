from dash import dcc, html
import dash_bootstrap_components as dbc
import dash_draggable
from guanaco.pages.single_cell.config import common_config
## Tab2: Violin2
def generate_violin_layout(adata, default_gene_markers,discrete_label_list,prefix):
    # Pre-create some common dropdowns and toggles
    violin1_transformation_selection = dbc.RadioItems(
        id=f'{prefix}-violin-log-or-zscore',
        options=[
            {'label': 'None', 'value': 'False'},
            {'label': 'Log', 'value': 'log'}
        ],
        value='log',
        inline=True,
        style={'marginLeft': '10px'}
    )

    violin2_transformation_selection = dbc.RadioItems(
        id=f'{prefix}-violin2-log-or-zscore',
        options=[
            {'label': 'None', 'value': 'False'},
            {'label': 'Log', 'value': 'log'}
        ],
        value='log',
        inline=False
    )

    violin_show_box1 = dbc.Checklist(
        id=f'{prefix}-show-box1',
        options=[{'label': 'Show Box Plot', 'value': 'show'}],
        value=[],
        switch=True
    )

    violin_show_scatter1 = dbc.Checklist(
        id=f'{prefix}-show-scatter1',
        options=[{'label': 'Show Scatter Points', 'value': 'show'}],
        value=[],
        switch=True
    )

    violin_show_box2 = dbc.Checklist(
        id=f'{prefix}-show-box2',
        options=[{'label': 'Show Box Plot', 'value': 'show'}],
        value=[],
        switch=True
    )

    violin_show_scatter2 = dbc.Checklist(
        id=f'{prefix}-show-scatter2',
        options=[{'label': 'Show Scatter Points', 'value': 'show'}],
        value=[],
        switch=True
    )

    # Meta1 dropdown (primary metadata)
    meta1_selection = dcc.Dropdown(
        id=f'{prefix}-meta1-selection',
        options=[{'label': meta, 'value': meta} for meta in discrete_label_list],
        value=discrete_label_list[0],
        clearable=False,
        placeholder="Select primary metadata"
    )

    # Meta2 dropdown (secondary metadata - optional)
    meta2_selection = dcc.Dropdown(
        id=f'{prefix}-meta2-selection',
        options=[{'label': 'None', 'value': 'none'}] + [{'label': meta, 'value': meta} for meta in discrete_label_list],
        value='none',
        clearable=False,
        placeholder="Select secondary metadata (optional)"
    )

    # Mode selection dropdown
    mode_selection = dcc.Dropdown(
        id=f'{prefix}-mode-selection',
        options=[
            {'label': 'Mode 1: One metadata only', 'value': 'mode1'},
            {'label': 'Mode 2: Facet by meta1, compare meta2', 'value': 'mode2'},
            {'label': 'Mode 3: Linear model (meta1 + confounder)', 'value': 'mode3'},
            {'label': 'Mode 4: Mixed model (meta1 + random effect)', 'value': 'mode4'},
        ],
        value='mode1',
        clearable=False
    )

    # Test method dropdown
    test_method_selection = dcc.Dropdown(
        id=f'{prefix}-test-method-selection',
        options=[
            {'label': 'None', 'value': 'none'},
            {'label': 'Mann-Whitney U', 'value': 'mwu-test'},
            {'label': 'T-test', 'value': 'ttest'},
            {'label': 'Kruskal-Wallis', 'value': 'kw-test'},
            {'label': 'ANOVA', 'value': 'anova'},
            {'label': 'Linear Model', 'value': 'linear-model'},
            {'label': 'Mixed Model', 'value': 'mixed-model'}
        ],
        value='none',
        clearable=False
    )

    violin2_gene_selection = dcc.Dropdown(
        id=f'{prefix}-violin2-gene-selection',
        options=[{'label': gene, 'value': gene} for gene in default_gene_markers],
        value=default_gene_markers[0] if default_gene_markers else None
    )

    group_selection = dcc.Dropdown(
        id=f'{prefix}-violin2-group-selection',
        multi=True,
        style={'font-size': '12px'}
    )

    # Final layout
    violin_layout = html.Div([
        # Violin Plot 1
        html.Div([
            html.Div([
                html.Label("Transformation:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                violin1_transformation_selection
            ], style={'marginBottom': '15px'}),
            html.Div([
                html.Label("More Options:", style={'fontWeight': 'bold'}),
                html.Div([violin_show_box1, violin_show_scatter1],
                        style={'display': 'flex', 'gap': '20px'})
            ], style={'marginBottom': '15px'}),
            # Add cache store for violin plot optimization
            dcc.Store(id=f'{prefix}-violin-plot-cache-store'),
            dcc.Loading(
                id="loading-violin1",
                type="circle",
                children=[
                    dcc.Graph(
                        id=f'{prefix}-violin-plot1',
                        config=common_config,
                        style={'width': '100%', 'minHeight': '400px'}
                    )
                ]
            ),
        ], style={'marginBottom': '30px', 'padding': '10px'}),

        # Violin Plot 2
        html.Div([
            html.H4("Split Violin/Grouped Violin", style={'textAlign': 'center', 'margin': '10px 0', 'fontWeight': 'bold'}),
            html.Label("Select Gene", style={'fontWeight': 'bold'}),
            violin2_gene_selection,
            html.Div([
                html.Label("Transformation:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                violin2_transformation_selection
            ], style={'marginBottom': '10px'}),
            html.Div([
                html.Div([
                    html.Label("Obs1 (Primary):", style={'fontWeight': 'bold'}),
                    meta1_selection
                ], style={'flex': '1'}),
                html.Div([
                    html.Label("Obs2 (Secondary):", style={'fontWeight': 'bold'}),
                    meta2_selection
                ], style={'flex': '1'}),
            ], style={'display': 'flex', 'marginBottom': '10px', 'gap': '10px'}),

            html.Div([
                html.Div([
                    html.Label("Analysis Mode:", style={'fontWeight': 'bold'}),
                    mode_selection
                ], style={'flex': '1'}),
                html.Div([
                    html.Label("Statistical Test:", style={'fontWeight': 'bold'}),
                    test_method_selection
                ], style={'flex': '1'}),
            ], style={'display': 'flex', 'marginBottom': '10px', 'gap': '10px'}),
            
            # Mode explanation helper text
            html.Div(
                id=f'{prefix}-mode-explanation',
                style={'fontSize': '12px', 'color': 'gray', 'marginBottom': '10px', 'fontStyle': 'italic'}
            ),

            html.Label("More Options:", style={'fontWeight': 'bold'}),
            violin_show_box2,
            violin_show_scatter2,

            # Draggable Graph for Violin 2
            dcc.Loading(
                id="loading-violin2",
                type="circle",
                children=[
                    html.Div([
                        dash_draggable.GridLayout(
                            id=f'{prefix}-draggable-violin2',
                            className='grid-layout-no-border',
                            children=[
                                html.Div(
                                    children=dcc.Graph(
                                        id=f'{prefix}-violin-plot2',
                                        config=common_config,
                                        responsive=True,
                                        style={
                                            "min-height": "0",
                                            "flex-grow": "1",
                                            "border": "none",
                                            "box-shadow": "none"
                                        }
                                    ),
                                    style={
                                        "height": '100%',
                                        "width": '100%',
                                        "display": "flex",
                                        "flex-direction": "column",
                                        "flex-grow": "0",
                                        "border": "none",
                                        "box-shadow": "none"
                                    }
                                )
                            ],
                            isResizable=True,
                            isDraggable=True,
                            height=30,
                            gridCols=12,
                            style={
                                "min-width": "300px",   # optional: prevents collapse
                                "min-height": "300px",  # optional: ensures initial size
                                "background": "transparent",
                                "padding": "0px"
                            }
                        )
                    ])
                ],
                style={
                    'padding': '10px',
                    'display': 'flex',
                    'flex-direction': 'column',
                    'height': '100%'
                }
            ),

                            
            # dash_draggable.GridLayout(
            #     id=f'{prefix}-draggable-violin2',
            #     className='grid-layout-no-border',
            #     children=[
            #         html.Div(
            #             dcc.Loading(
            #                 id="loading-violin2",
            #                 type="circle",
            #                 children=[
            #                     dcc.Graph(
            #                         id=f'{prefix}-violin-plot2',
            #                         config=common_config,
            #                         responsive=True,
            #                         style={"minHeight": "300px", "flex-grow": "1"}
            #                     )
            #                 ],
            #                 style={
            #                     "marginTop": "10px",
            #                     "height": "100%",
            #                     "display": "flex",
            #                     "flex-direction": "column"
            #                 }
            #             ),
            #             style={
            #                 "height": '100%',
            #                 "width": '100%',
            #                 "display": "flex",
            #                 "flex-direction": "column"
            #             }
            #         ),
            #     ],
            #     isResizable=True,
            #     isDraggable=True,
            #     height=30,
            #     ncols=12,
            #     style={'background': 'transparent', 'padding': '0px', 'minHeight': '400px', 'flexGrow': '0'}
            # ),

            html.Label("Labels in Group:", style={'fontWeight': 'bold', 'marginTop': '10px'}),
            group_selection,
        ], style={'padding': '10px'})
    ], style={'width': '100%'})



    return violin_layout

