import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
from guanaco.data_loader import color_config

def plot_stacked_bar(x_meta, y_meta, norm, adata, color_map=None, y_order=None):
    """
    Plot stacked bar chart.
    
    Arguments:
        x_meta: Column name for x-axis groups
        y_meta: Column name for stacked layers (color)
        norm: 'prop' for proportions, 'count' for counts
        adata: AnnData object
        color_map: Color mapping dictionary or list
        y_order: List specifying order of stacked layers (from Select Labels)
    """
    # Create count dataframe
    count_df = adata.obs.groupby([x_meta, y_meta]).size().reset_index(name='count')
    count_df[x_meta] = count_df[x_meta].astype(str)
    count_df[y_meta] = count_df[y_meta].astype(str)

    # Calculate proportions if needed
    if norm == 'prop':
        count_df['prop'] = count_df.groupby(x_meta)['count'].transform(lambda x: x / x.sum())
        y_value = 'prop'
        y_label = 'Proportion'
    else:
        y_value = 'count'
        y_label = 'Cell Count'
    
    # Set up color mapping - use the provided color_map to ensure consistency
    if isinstance(color_map, dict):
        color_discrete_map = color_map
    elif color_map is None:
        categories = sorted(count_df[y_meta].unique())
        predefined_colors = color_config
        color_discrete_map = {cat: predefined_colors[i % len(predefined_colors)] for i, cat in enumerate(categories)}
    else:
        categories = sorted(count_df[y_meta].unique())
        color_discrete_map = {cat: color_map[i % len(color_map)] for i, cat in enumerate(categories)}

    # Apply y_order if specified (this comes from Select Labels)
    if y_order is not None and len(y_order) > 0:
        # Convert y_order to strings for consistency
        y_order_str = [str(y) for y in y_order]
        # Filter the dataframe to only include categories in y_order
        count_df = count_df[count_df[y_meta].isin(y_order_str)]
        # Sort the dataframe to match the order in y_order
        count_df[y_meta] = pd.Categorical(count_df[y_meta], categories=y_order_str, ordered=True)
        count_df = count_df.sort_values([x_meta, y_meta])
        category_orders = {y_meta: y_order_str}
    else:
        # If no order specified, use all categories in sorted order
        category_orders = None

    # Create the plot
    fig = px.bar(
        count_df,
        x=x_meta,
        y=y_value,
        color=y_meta,
        labels={x_meta: f'{x_meta}', y_value: y_label, y_meta: f'{y_meta}'},
        barmode='stack',
        color_discrete_map=color_discrete_map,
        category_orders=category_orders  # This ensures the stacking order follows y_order
    )

    # Update layout
    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        xaxis=dict(
            showgrid=False, 
            title_font=dict(size=18),
            categoryorder='total ascending'
        ),
        yaxis=dict(showgrid=False, title_font=dict(size=18)),
        legend=dict(
            title=y_meta,
            orientation="v",
            yanchor="middle",
            y=0.5,
            xanchor="left",
            x=1.02
        ),
        margin=dict(r=150)  # Make room for legend
    )

    fig.update_xaxes(showline=True, linewidth=2, linecolor='black')
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black')

    return fig