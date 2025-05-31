from dash import dcc, html
import dash_bootstrap_components as dbc
import dash_draggable
from guanaco.config import common_config

def generate_stacked_bar_layout(discrete_label_list,prefix):

    x_meta_dropdown = html.Div([
        html.Label("Cell info 1(X-axis):", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
        dcc.Dropdown(
            id=f'{prefix}-x-meta',
            options=[{'label': meta, 'value': meta} for meta in discrete_label_list],
            value=discrete_label_list[0],
            clearable=False,
            style={'marginBottom': '15px'}
        )
    ], style={'flex': '1'})

    y_meta_dropdown = html.Div([
        html.Label("Cell info 2(colour by):", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
        dcc.Dropdown(
            id=f'{prefix}-y-meta',
            options=[{'label': meta, 'value': meta} for meta in discrete_label_list],
            value=discrete_label_list[1] if len(discrete_label_list) > 1 else 'none',
            clearable=False,
            style={'marginBottom': '15px'}
        )
    ], style={'flex': '1'})

    norm_box = html.Div([
        html.Label("Plot value:", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
        dbc.RadioItems(
            id=f'{prefix}-norm-box',
            options=[
                {'label': 'Proportion', 'value': 'prop'},
                {'label': 'Count', 'value': 'count'},
            ],
            value='prop',
            inline=True
        )
    ], {'padding': '20px'})
    draggable_bar = html.Div([
    dash_draggable.GridLayout(
        id=f'{prefix}-draggable',
        className='grid-layout-no-border',
        children=[
            html.Div(children=[
                dcc.Graph(
                    id=f'{prefix}-stacked-bar-plot',
                    config=common_config,
                    responsive=True,
                    style={
                        "min-height": "0",
                        "flex-grow": "1",
                    }
                )
            ],
            style={
                "height": '100%',
                "width": '100%',
                "display": "flex",
                "flex-direction": "column",
                "flex-grow": "0",
            }),
        ],
        isResizable=True,
        isDraggable=True,
        height=30,
        gridCols=12,
    )
])


    bar_layout = html.Div([
        html.Div(
            [x_meta_dropdown, y_meta_dropdown],
            style={'display': 'flex', 'justify-content': 'space-between'}
        ),
        norm_box,
        draggable_bar
    ], style={'padding': '20px',
               'marginBottom': '15px',
               })

    return bar_layout
