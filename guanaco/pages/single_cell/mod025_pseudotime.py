from dash import dcc, html
import dash_bootstrap_components as dbc


def generate_pseudotime_layout(prefix):
    """Generate the layout for the pseudotime plot tab."""
    
    layout = html.Div([
        # Controls container
        dbc.Row([
            dbc.Col([
                # Min expression threshold slider
                html.Div([
                    html.Label('Minimum Expression Threshold:', 
                              className="control-label"),
                    dcc.Slider(
                        id=f'{prefix}-pseudotime-min-expr-slider',
                        min=0,
                        max=5,
                        step=0.1,
                        value=0.5,
                        marks=None,
                        tooltip={"placement": "bottom", "always_visible": True},
                        className="dbc-slider"
                    )
                ], style={'marginBottom': '20px'}),
                
                # Transformation selection
                html.Div([
                    html.Label('Transformation:', 
                              style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                    dbc.RadioItems(
                        id=f'{prefix}-pseudotime-transformation',
                        options=[
                            {'label': 'None', 'value': 'none'},
                            {'label': 'Log', 'value': 'log'},
                            {'label': 'Z-score', 'value': 'z_score'}
                        ],
                        value='none',
                        inline=True,
                        style={'fontSize': '14px'}
                    )
                ], style={'marginBottom': '20px'}),
                
                # Pseudotime key selection (optional, in case there are multiple)
                html.Div([
                    html.Label('Pseudotime Column:', 
                              style={'fontWeight': 'bold', 'marginBottom': '5px'}),
                    dcc.Dropdown(
                        id=f'{prefix}-pseudotime-key-dropdown',
                        options=[],  # Will be populated dynamically
                        value=None,  # Will be set dynamically
                        placeholder='Select pseudotime column',
                        style={'marginBottom': '15px'},
                        clearable=False
                    )
                ], style={'marginBottom': '20px'}),
                
            ], width=12)
        ], style={'marginBottom': '20px', 'padding': '15px', 'backgroundColor': '#f8f9fa', 'borderRadius': '5px'}),
        
        # Plot container
        dbc.Row([
            dbc.Col([
                dcc.Loading(
                    id=f"{prefix}-pseudotime-loading",
                    type="default",
                    children=[
                        dcc.Graph(
                            id=f'{prefix}-pseudotime-plot',
                            style={'height': 'auto', 'minHeight': '600px'},
                            config={
                                'displayModeBar': True,
                                'displaylogo': False,
                                'modeBarButtonsToRemove': ['lasso2d', 'select2d'],
                                'toImageButtonOptions': {
                                    'format': 'png',
                                    'filename': 'pseudotime_plot',
                                    'height': 800,
                                    'width': 1200,
                                    'scale': 2
                                }
                            }
                        )
                    ]
                )
            ], width=12)
        ])
    ])
    
    return layout