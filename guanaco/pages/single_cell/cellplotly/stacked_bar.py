import plotly.express as px
from guanaco.data_loader import color_config

def plot_stacked_bar(x_meta, y_meta, norm, adata, color_map=None):
    count_df = adata.obs.groupby([x_meta, y_meta]).size().reset_index(name='count')
    count_df[x_meta] = count_df[x_meta].astype(str)
    count_df = count_df.sort_values(by=x_meta)

    if norm == 'prop':
        count_df['prop'] = count_df.groupby(x_meta)['count'].transform(lambda x: x / x.sum())
        y_value = 'prop'
        y_label = 'Proportion'
    else:
        y_value = 'count'
        y_label = 'Cell Count'
    
    # Use provided color_map directly if it's a dictionary (fixed mapping)
    # Otherwise create mapping based on current categories
    if isinstance(color_map, dict):
        # Use the fixed color mapping passed from the callback
        color_discrete_map = color_map
    elif color_map is None:
        # Fallback to default colors for current categories
        categories = sorted(count_df[y_meta].unique())
        predefined_colors = color_config
        color_discrete_map = {cat: predefined_colors[i % len(predefined_colors)] for i, cat in enumerate(categories)}
    else:
        # color_map is a list of colors
        categories = sorted(count_df[y_meta].unique())
        color_discrete_map = {cat: color_map[i % len(color_map)] for i, cat in enumerate(categories)}

    # 4. Create plot with the color map
    fig = px.bar(
        count_df,
        x=x_meta,
        y=y_value,
        color=y_meta,
        labels={x_meta: f'{x_meta}', y_value: y_label, y_meta: f'{y_meta}'},
        barmode='stack',
        color_discrete_map=color_discrete_map
    )

    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        xaxis=dict(showgrid=False, title_font=dict(size=18)),
        yaxis=dict(showgrid=False, title_font=dict(size=18)),
    )

    fig.update_xaxes(showline=True, linewidth=2, linecolor='black')
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black')

    return fig
