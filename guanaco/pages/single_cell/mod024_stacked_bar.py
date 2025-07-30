from dash import dcc, html
import dash_bootstrap_components as dbc
import dash_draggable
import dash_ag_grid as dag
from guanaco.config import common_config

def generate_stacked_bar_layout(discrete_label_list, prefix):
    """
    Generate layout for stacked bar plot.
    
    Arguments:
        discrete_label_list: List of discrete metadata columns
        prefix: Unique ID prefix for this instance
    """
    
    # Dropdown to select x-axis metadata
    x_axis_dropdown = html.Div([
        html.Label("Cell info (x-axis):", style={'fontWeight': 'bold', 'marginBottom': '5px'}),
        dcc.Dropdown(
            id=f'{prefix}-stacked-bar-x-axis',
            options=[{'label': meta, 'value': meta} for meta in discrete_label_list],
            value=discrete_label_list[0] if discrete_label_list else None,
            clearable=False,
            placeholder="Select metadata for x-axis",
            style={'marginBottom': '15px'}
        )
    ], style={'marginBottom': '15px'})
    
    # Normalization option
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
    ], style={'marginBottom': '15px'})
    
    # Draggable stacked bar plot
    draggable_bar = dash_draggable.GridLayout(
        id=f'{prefix}-draggable',
        className='grid-layout-no-border',
        children=[
            html.Div(children=[
                dcc.Graph(
                    id=f'{prefix}-stacked-bar-plot',
                    config=common_config,
                    responsive=True,
                    style={
                        "min-height": "400px",
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
    
    # Info text
    info_text = html.Div([
        html.P([
            html.I(className="fas fa-info-circle", style={'marginRight': '5px'}),
            "The stacked layers come from 'Select Annotation' and 'Select Labels' in the left control panel"
        ], style={'color': '#6c757d', 'fontSize': '14px', 'marginBottom': '15px'})
    ])
    
    # Draggable AG Grid component for x-axis ordering (will be placed after plot)
    x_axis_order_component = html.Div([
        html.Label("X-axis group order:", 
                   style={'fontWeight': 'bold', 'marginBottom': '10px'}),
        dag.AgGrid(
            id=f'{prefix}-stacked-bar-x-order-grid',
            rowData=[],  # No rows, only headers
            columnDefs=[],  # Will be populated dynamically based on x-axis selection
            defaultColDef={
                "sortable": False,
                "filter": False,
                "resizable": True,
                "suppressMenu": True,
                "headerHeight": 40,
                "minWidth": 120,
                "width": 150,
                "headerClass": "ag-header-cell-center"
            },
            dashGridOptions={
                "headerHeight": 40,
                "rowHeight": 0,
                "suppressRowClickSelection": True,
                "suppressCellSelection": True,
                "suppressMovableColumns": False,  # Enable column dragging
                "animateRows": False,
                "suppressHorizontalScroll": False,
                "onColumnMoved": True,  # Enable column moved event
                "suppressLoadingOverlay": True,
                "suppressNoRowsOverlay": True,
                "suppressDisplayTotal": True
            },
            style={"height": "40px", "marginBottom": "10px", "overflow": "hidden"},
            className="ag-theme-alpine"
        ),
        html.P("Drag column headers to reorder x-axis groups.", 
               style={'fontSize': '12px', 'color': '#6c757d', 'marginTop': '5px', 'marginBottom': '15px'})
    ])
    
    # Store to track column order
    column_order_store = dcc.Store(id=f'{prefix}-x-axis-column-order-store', data=[])
    
    # Final layout
    bar_layout = html.Div([
        column_order_store,
        info_text,
        x_axis_dropdown,
        norm_box,
        draggable_bar,
        x_axis_order_component  # Moved to bottom after the plot
    ], style={
        'padding': '20px',
        'marginBottom': '15px',
    })

    return bar_layout