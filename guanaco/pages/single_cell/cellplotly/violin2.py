import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from scipy.stats import ttest_ind, mannwhitneyu, f_oneway, kruskal
import statsmodels.api as sm
from statsmodels.formula.api import ols
from guanaco.data_loader import color_config
def calculate_p_values(df, method, groupby, hue,labels):
    if method == 'None':
        return None

    p_values = {}
    unique_hues = df[hue].unique()

    if method in ['kw-test', 'anova']:
        for condition in unique_hues:
            condition_data = df[df[hue] == condition]
            groups = [group['Expression'].values for _, group in condition_data.groupby(groupby, observed=True)]
            test_func = kruskal if method == 'kw-test' else f_oneway
            _, p_val = test_func(*groups)
            p_values[condition] = p_val

    elif method in ['mwu-test', 'ttest']:
        for cell_type in labels:
            cell_data = df[df[groupby] == cell_type]
            group_1 = cell_data[cell_data[hue] == unique_hues[0]]['Expression']
            group_2 = cell_data[cell_data[hue] == unique_hues[1]]['Expression']
            test_func = mannwhitneyu if method == 'mwu-test' else ttest_ind
            kwargs = {'alternative': 'two-sided'} if method == 'mwu-test' else {'equal_var': False}
            _, p_val = test_func(group_1, group_2, **kwargs)
            p_values[cell_type] = p_val

    elif method == 'two-way-anova':
        formula = f'Expression ~ C({groupby}) * C({hue})'
        model = ols(formula, data=df).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)
        interaction_term = f'C({groupby}):C({hue})'
        p_values['interaction'] = anova_table.loc[interaction_term, 'PR(>F)']

    return p_values


def assign_colors(hue_levels, hue_color_map=None):
    if hue_color_map is None:
        hue_color_map = {}
    default_colors = px.colors.qualitative.Plotly
    for i, level in enumerate(sorted(hue_levels)):
        if level not in hue_color_map:
            hue_color_map[level] = color_config[i] if i < len(color_config) else default_colors[i % len(default_colors)]
    return hue_color_map


def add_p_value_annotations(fig, p_values, method, df,labels):
    if not p_values:
        return

    if method in ['mwu-test', 'ttest']:
        for group, p_val in p_values.items():
            if labels is None or group in labels:
                fig.add_annotation(
                    x=group,
                    y=max(df['Expression']) * 1.05,
                    # text=f"{p_val:.3f}",
                    text = f"<b>{p_val:.3f}</b>" if p_val < 0.05 else f"{p_val:.3f}",
                    showarrow=False,
                    font=dict(size=10, color='DarkSlateGrey')
                )
    else:
        p_text = "\n".join([f"{k}: <b>{v:.3g}</b>" for k, v in p_values.items()])
        fig.add_shape(
            type="path",
            path="M 0.02 ,1.05 L 0.02,1.09 L 0.98,1.09 L 0.98,1.05",
            xref="paper",
            yref="paper",
            line=dict(color="DarkSlateGrey", width=2)
        )
        fig.add_annotation(
            text=p_text,
            align="center",
            showarrow=False,
            xref="paper",
            yref="paper",
            x=0.5,
            y=1.20,
            font=dict(size=10, color="DarkSlateGrey"),
            bgcolor="rgba(255, 255, 255, 0)"
        )


def plot_violin2(adata, key, groupby, hue, transformation,
                 show_box=True, show_points=True, p_value=None,labels=None,
                  hue_color_map=None):

    if not hue or not groupby or not key:
        return go.Figure()


    # Prepare data
    if key not in adata.var_names:
        return go.Figure()
    
    if key not in adata.var_names:
        return go.Figure()

    expression_data = adata[:, key].X.toarray().flatten()
    if transformation == 'log':
        expression_data = np.log1p(expression_data)
    elif transformation == 'z_score':
        expression_data = (expression_data - expression_data.mean()) / expression_data.std()

    df = pd.DataFrame({
        'Expression': expression_data,
        groupby: adata.obs[groupby],
        hue: adata.obs[hue],
        'Feature': key
    })

    hue_levels = df[hue].unique()
    hue_color_map = assign_colors(hue_levels, hue_color_map)
    points_mode = 'all' if show_points else False

    fig = go.Figure()
    grouped_data = df.groupby(groupby,observed=True)
    if len(hue_levels) == 2:  # Split violin
        for label in labels:
            group_data = grouped_data.get_group(label)
            for hue_val in sorted(hue_levels):
                side = 'negative' if hue_val == hue_levels[0] else 'positive'
                sub_df = group_data[group_data[hue] == hue_val]
                fig.add_trace(go.Violin(
                    x=sub_df[groupby],
                    y=sub_df['Expression'],
                    legendgroup=hue_val,
                    width=0.8,
                    scalemode='width',  # uniform width
                    scalegroup=str(label),
                    name=hue_val,
                    box_visible=show_box,
                    points=points_mode,
                    bandwidth=0.2,
                    spanmode='hard',
                    meanline_visible=True,
                    side=side,
                    line_color='DarkSlateGrey',
                    fillcolor=hue_color_map[hue_val],
                    hoveron='violins',
                    hoverinfo='y',
                    showlegend=(label == labels[0]),
                    jitter=0.05,
                ))
        p_values = calculate_p_values(df, p_value, groupby, hue, labels=labels)
        add_p_value_annotations(fig, p_values, p_value, df, labels=labels)
        fig.update_layout(violinmode='overlay')

    else:  # grouped violin
        for hue_val in sorted(hue_levels):
            subset = df[df[hue] == hue_val]
            fig.add_trace(go.Violin(
                y=subset['Expression'],
                x=subset[groupby],
                # width=0.8,
                name=hue_val,
                meanline_visible=True,
                points=points_mode,
                box_visible=show_box,
                bandwidth=0.2,
                span=[0, max(df['Expression'])],
                legendgroup=hue_val,
                fillcolor=hue_color_map[hue_val],
                line_color='DarkSlateGrey',
                hoveron='violins',
                hoverinfo='y',
                spanmode='hard',
                jitter=0.05,
            ))
            fig.update_layout(violingap=0.15,violinmode='group')

    # Shared layout updates
    fig.update_layout(
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(size=12, color='DarkSlateGrey'),
        title=dict(
            text=f"Expression of {key} Across {groupby}",
            x=0.5,
            xanchor='center',
            y=0.9,
            font=dict(size=16, color='DarkSlateGrey')
        ),
        yaxis_range=[0, max(df['Expression']) * 1.05],
        xaxis_title=groupby,
        yaxis_title='Expression',
            legend=dict(
        font=dict(size=10) )
    )

    if show_points:
        fig.update_traces(marker=dict(size=1.5), selector=dict(type='violin'))

    fig.update_yaxes(
        showgrid=False,
        showline=True, linewidth=2, linecolor='black',
        title=dict(text=f"{key}", font=dict(size=14, color='DarkSlateGrey'))
    )
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black')

    return fig