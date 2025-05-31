from dash import dcc, html
import dash_bootstrap_components as dbc
import dash_draggable
from guanaco.config import common_config
## Tab2: Violin2
def generate_violin_layout(adata, default_gene_markers,discrete_label_list,prefix):
    # Pre-create some common dropdowns and toggles
    violin1_transformation_selection = dbc.RadioItems(
        id=f'{prefix}-violin-log-or-zscore',
        options=[
            {'label': 'None', 'value': 'False'},
            {'label': 'Log', 'value': 'log'},
            {'label': 'Z-score(across cell)', 'value': 'z_score'}
        ],
        value='log',
        inline=True,
        style={'marginLeft': '10px'}
    )

    violin2_transformation_selection = dbc.RadioItems(
        id=f'{prefix}-violin2-log-or-zscore',
        options=[
            {'label': 'None', 'value': 'False'},
            {'label': 'Log', 'value': 'log'},
            {'label': 'Z-score', 'value': 'z_score'}
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

    comparison_selection = dcc.Dropdown(
        id=f'{prefix}-p-value-selection',
        options=[
            {'label': 'None', 'value': 'None'},
            {'label': 'Within Group', 'value': 'binary'},
            {'label': 'Between Group', 'value': 'multi'},
        ],
        value='binary',
        clearable=False
    )

    binary_selection = dcc.Dropdown(
        id=f'{prefix}-binary-selection',
        options=[{'label': meta, 'value': meta} for meta in discrete_label_list],
        value=discrete_label_list[0],
        clearable=False
    )

    multi_class_selection = dcc.Dropdown(
        id=f'{prefix}-multi-class-selection',
        options=[{'label': meta, 'value': meta} for meta in discrete_label_list],
        value=discrete_label_list[len(discrete_label_list)//2],
        clearable=False
    )

    p_value_method = dcc.Dropdown(
        id=f'{prefix}-p-value-method',
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
                    html.Label("Group by:", style={'fontWeight': 'bold'}),
                    multi_class_selection
                ], style={'flex': '1'}),
                html.Div([
                    html.Label("Splitted by:", style={'fontWeight': 'bold'}),
                    binary_selection
                ], style={'flex': '1'}),
            ], style={'display': 'flex', 'marginBottom': '10px', 'gap': '10px'}),

            html.Div([
                html.Div([
                    html.Label("Comparison Type:", style={'fontWeight': 'bold'}),
                    comparison_selection
                ], style={'flex': '1'}),
                html.Div([
                    html.Label("Test Method:", style={'fontWeight': 'bold'}),
                    p_value_method
                ], style={'flex': '1'}),
            ], style={'display': 'flex', 'marginBottom': '10px', 'gap': '10px'}),

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

