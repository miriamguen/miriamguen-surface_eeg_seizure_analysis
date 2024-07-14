"This "
# imports
import json
import os
from os.path import join
from pathlib import Path

import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scikit_posthocs as sp
import scipy.stats as ss
import seaborn as sns
import statsmodels.api as sa
import statsmodels.formula.api as sfa
from sklearn.feature_selection import f_classif

# define data version and paths
DATA_ALIAS = 'example_data'
# folder definition
BASE_PATH = os.getcwd().split('feature_analysis')[0]
DATA_DIR = os.path.join(BASE_PATH, DATA_ALIAS, 'features')
FIG_DIR = os.path.join(DATA_DIR, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

onset_palette_only_bif = {'SN/SubH': "#7491fb", 'SupH': "#a9d39e", 'SNIC': "#fdfd96"}
onset_palette_seizure = {"Unclear": 'navy', **onset_palette_only_bif}
onset_palette = {'Control': '#4fd6ab', **onset_palette_seizure}

offset_palette_only_bif = {'FLC': "#89e4db", 'SH/SNIC': "#fbb374", 'SupH': "#a9d39e"}
offset_palette_seizure = {"Unclear": 'navy', **offset_palette_only_bif}
offset_palette = {'Control': '#4fd6ab', **offset_palette_seizure}

group_palette = {"Clear": "#8ad9e3", "Unclear": "navy", 'Control': '#4fd6ab'}
palette_clear = {'Clear': "#8ad9e3", 'Unclear': 'navy'}
cmap = 'winter'  # 'jet'#'inferno'#'Spectral'#'cool'#'viridis'#'Spectral'#'PRGn'#, 'inferno'
plt.style.reload_library()
plt.style.use(['science'])
params = {"font.family": 'sans-serif',
          "font.sans-serif": 'calibri',
          'text.usetex': False,
          'axes.titlesize': 40,
          'axes.labelsize': 28,
          'axes.grid': True,
          'axes.grid.axis': 'y',
          'axes.grid.which': 'major',
          'axes.edgecolor': 'black',
          'axes.facecolor': 'white',
          'axes.linewidth': 0.4,
          'ytick.labelsize': 22,
          'xtick.labelsize': 22,
          'legend.facecolor': 'white',
          'legend.fontsize': 18,
          'figure.facecolor': 'white',
          'savefig.transparent': True,
          'backend': 'svg'
          }
pylab.rcParams.update(params)
sns.set(rc=dict(pylab.rcParams))


# help functions data cleaning and feature selection
def arrange_tables(feature_table, common_components, start_end):
    labels = common_components[common_components[f'clear_{start_end}'] == 1]
    labels = labels.replace('not_clear', 'Unclear')
    table = pd.concat([feature_table.set_index(['file_name', 'comp_num'], drop=True),
                       labels.set_index(['file_name', 'comp_num'], drop=True)], join='outer', axis=1)

    table['clear_start'] = table[f'clear_start'].apply(lambda x: 'Clear' if x == 1 else "Unclear")
    table['clear_end'] = table[f'clear_end'].apply(lambda x: 'Clear' if x == 1 else "Unclear")
    table.bif_start.fillna("Unclear", inplace=True)
    table.bif_end.fillna("Unclear", inplace=True)
    table.scm.fillna(False, inplace=True)
    table.clear_start.fillna(0, inplace=True)
    table.clear_end.fillna(0, inplace=True)
    table.fillna(0, inplace=True)
    table.replace([np.inf], 0, inplace=True)

    table['max_score'] = table[['brain', 'mara']].max(axis=1)
    # mark the change points detected
    table.loc[:, 'Change point'] = table['bif_time'].apply(lambda x: 'Detected' if x > 0 else 'Undetected')
    table['change_point'] = table['Change point'].apply(lambda x: 1 if x == 'Detected' else 0)
    return table


def aggregate_control(data, control_data, bif_time=-np.inf, score_th=0):
    control_data = control_data[control_data['bif_time'] > bif_time]
    control_data = control_data[control_data['score'] > score_th]
    # control_data = control_data[control_data['rms_error']<5e-6]
    control_data['clear_start'] = 'Control'
    control_data['clear_end'] = 'Control'

    data = data[data['bif_time'] > bif_time]
    data = data[data['score'] > score_th]
    # control_data = control_data[control_data['rms_error']<5e-6]

    plot_data = pd.concat([data.reset_index(), control_data.reset_index()])

    return plot_data.reset_index(drop=True)


def plot_boxes(feature, table, group, ax):
    bifs = sorted(list(set(table[group])), key=str.lower)

    data = [table.loc[table[group] == bif, feature] for bif in bifs]
    n = [len(x) for x in data]
    ax.set_title(f' {feature.capitalize()} by {group}')
    # Creating axes instance

    # Creating plot
    bp = ax.boxplot(data, notch='True')
    # changing color and line width of medians
    for median in bp['medians']:
        median.set(color='teal',
                   linewidth=4)

    # changing color and line width of caps
    for cap in bp['caps']:
        cap.set(color='#8B008B',
                linewidth=3)

    # changing color and line width of whiskers
    for whisker in bp['whiskers']:
        whisker.set(color='#8B008B',
                    linewidth=2,
                    linestyle=":")

    # x-axis labels
    formatted_ticks = [f'{b}\n n={len(x)}' for b, x in zip(bifs, data)]
    ax.set_xticklabels(formatted_ticks, fontsize=12, fontweight='bold')
    return ax


def plot_scatter(data, feature1, feature2, feature3, palette='coolwarm', alpha=0.8, sizes=1):
    figure, ax = plt.subplots(1, figsize=(10, 10))
    ax = sns.scatterplot(x=feature1, y=feature2, hue=feature3, data=data, palette=palette, ax=ax, legend='auto',
                         alpha=alpha, sizes=sizes)
    ax.set_title(f'{feature1} vs {feature2} vs {feature3}: all')

    figure.tight_layout(pad=0.5)
    plt.show()
    return figure, ax


def plot_violin(feature, table, group, palette, ax):
    # order = sorted(list(set(palette.keys())),key=str.lower)
    ax = sns.violinplot(x=group, y=feature, data=table, palette=palette, order=list(palette.keys()), ax=ax)
    ax.set_title(f" {feature.replace('_', ' ').capitalize()} by {group.replace('_', ' ')}");
    formatted_ticks = [f"{b.get_text()}\n n={sum(table[group] == b.get_text())}" for b in ax.get_xticklabels()]
    ax.set_xticklabels(formatted_ticks)
    ax.set_ylabel(feature)
    return ax


def plot_bar(feature, table, group, palette, ax, error_bars):
    edge_color = 0.5
    # order = sorted(list(set(palette.keys())),key=str.lower)
    # figure = plt.figure(figsize=(10,7))
    if error_bars:
        sns.barplot(x=group, y=feature, data=table, palette=palette, order=list(palette.keys()),
                    ax=ax, n_boot=5000)
    else:
        sns.barplot(x=group, y=feature, data=table, palette=palette, order=list(palette.keys()),
                    ax=ax, errorbar=None)

    ax.set_title(
        f" {feature.replace('_', ' ').capitalize()} by {group.replace('_', ' ').replace('bif', 'bifurcation:')}");
    formatted_ticks = [f"{b.get_text()}\n n={sum(table[group] == b.get_text())}" for b in ax.get_xticklabels()]
    ax.set_xticklabels(formatted_ticks)

    ax.set_ylabel(feature.capitalize())
    ax.set_xlabel(group.replace('_', ' ').replace('bif', 'bifurcation').capitalize())
    ax.grid(False)
    return ax


def plot_selected_one(data, feature, by, ax, title=None, kind='bar', palette=onset_palette, error_bars=True):
    if kind == 'bar':
        ax = plot_bar(feature, data, by, palette, ax, error_bars)
    else:
        ax = plot_violin(feature, data, by, palette, ax)

    if title == None:
        ax.set_title(f"{feature.replace('_', ' ')}")  # " p={anova.loc[feature, 'P']:.2f}")
    else:
        ax.set_title(f"{title} ")  # p={anova.loc[feature, 'P']:.2f}")

    ax.set_xlabel('')
    ax.set_ylabel(feature.replace('_', ' '))

    return ax


def plot_selected(data, features, anova, term, by, kind='bar', palette=onset_palette, error_bars=True):
    features_to_use = list(filter(lambda feature: term in feature, features))
    figure, axes = plt.subplots(len(features_to_use), 1, figsize=(7, 5 * len(features_to_use)))
    for i, feature in enumerate(features_to_use):
        if kind == 'bar':
            axes[i] = plot_bar(feature, data, by, palette, axes[i], error_bars)
        else:
            axes[i] = plot_violin(feature, data, by, palette, axes[i])
        axes[i].set_title(f"{feature.replace('_', ' ')} p = {anova.loc[feature, 'P']:.4f}", fontsize=20)
        axes[i].set_xlabel('Component Group')
        axes[i].set_ylabel(feature.replace('_', ' '))
    figure.tight_layout(pad=1)
    return figure, axes


def plot_feature_and_save(feature, data_onset, data_offset, onset_palette, offset_palette, FIG_DIR, ylim=None,
                          bif_clear='bif', error_bars=True):
    fig, ax = plt.subplots(1, 2, sharey='row', figsize=(20, 10))

    plot_selected_one(data_onset, feature, f'{bif_clear}_start', ax[0],
                      title='Onset ', kind='bar', palette=onset_palette,
                      error_bars=error_bars)

    plot_selected_one(data_offset, feature, f'{bif_clear}_end', ax[1],
                      title='Offset ', kind='bar', palette=offset_palette,
                      error_bars=error_bars)

    ax[0].set_ylabel(feature.replace('_', ' '))
    ax[1].set_ylabel('')
    if ~(ylim is None):
        ax[0].set_ylim(ylim)

    fig.suptitle(feature.replace('_', ' ').capitalize(), fontsize=40)
    fig.tight_layout(pad=1)

    fig_name = f'{feature}'
    print(fig_name)
    os.makedirs(FIG_DIR, exist_ok=True)
    fig.savefig(os.path.join(FIG_DIR, f'{fig_name}.svg'), transparent=True, bbox_inches='tight', pad_inches=1)


# Explore onset and offset features

# Load seizure data
file_list = os.listdir(join(DATA_DIR, 'seizure'))
onset_feature_files = list(filter(lambda x: '_start_' in x, file_list))
offset_feature_files = list(filter(lambda x: '_end_' in x, file_list))
overall_score_files = list(filter(lambda x: '_score_' in x, file_list))
len(overall_score_files)

seizure_score_table = pd.DataFrame()
if len(overall_score_files) > 0:
    for file_name in overall_score_files:
        score_file = pd.read_csv(join(DATA_DIR, 'seizure', file_name)).fillna(0)
        seizure_score_table = pd.concat([seizure_score_table, score_file], axis=0)

    seizure_score_table = seizure_score_table.set_index(['file_name', 'comp_num']).rename(lambda x: f'{x}_both', axis=1)

onset_feature_table = pd.DataFrame()

for file_name in onset_feature_files:
    onset_file = pd.read_csv(join(DATA_DIR, 'seizure', file_name)).fillna(0)
    onset_feature_table = pd.concat([onset_feature_table, onset_file], axis=0)

total_number = len(onset_feature_table)
success_number = len(onset_feature_table)
print(f'successes reading {success_number} out of {total_number} at onset')

offset_feature_table = pd.DataFrame()

for file_name in offset_feature_files:
    offset_file = pd.read_csv(join(DATA_DIR, 'seizure', file_name)).fillna(0)
    offset_feature_table = pd.concat([offset_feature_table, offset_file], axis=0)
total_number = len(offset_feature_table)
success_number = len(offset_feature_table)
print(f'successes reading {success_number} out of {total_number} at offset')

# load control data
file_list = os.listdir(join(DATA_DIR, 'control'))

onset_feature_files = list(filter(lambda x: '_start_' in x, file_list))
len(onset_feature_files)
offset_feature_files = list(filter(lambda x: '_end_' in x, file_list))
len(offset_feature_files)
overall_score_files = list(filter(lambda x: '_score_' in x, file_list))
len(overall_score_files)

control_score_table = pd.DataFrame()
if len(overall_score_files) > 0:
    for file_name in overall_score_files:
        score_file = pd.read_csv(join(DATA_DIR, 'control', file_name)).fillna(0)
        control_score_table = pd.concat([control_score_table, score_file], axis=0)

    control_score_table = control_score_table.set_index(['file_name', 'comp_num']).rename(lambda x: f'{x}_both', axis=1)

onset_control_feature_table = pd.DataFrame()
for file_name in onset_feature_files:
    onset_file = pd.read_csv(join(DATA_DIR, 'control', file_name)).fillna(0)
    onset_control_feature_table = pd.concat([onset_control_feature_table, onset_file], axis=0)

offset_control_feature_table = pd.DataFrame()
for file_name in offset_feature_files:
    offset_file = pd.read_csv(join(DATA_DIR, 'control', file_name)).fillna(0)
    offset_control_feature_table = pd.concat([offset_control_feature_table, offset_file], axis=0)

# Read a list with the matching components with unclear in case of mismatch
common_components = pd.read_csv(join(BASE_PATH, 'data', 'labels', 'revisions', 'common_components.csv'),
                                index_col='index')

common_components = common_components.drop(columns=['bif_start_r', 'bif_end_r', 'freq_r',
                                                    'scm_r', 'clear_start_r', 'clear_end_r', 'notes_r', 'freq_m'])

common_components = common_components.rename({x: x.split('_m')[0] for x in common_components.columns.values}, axis=1)

# add metadata to features
onset_seizure_table = arrange_tables(onset_feature_table, common_components, 'start')
offset_seizure_table = arrange_tables(offset_feature_table, common_components, 'end')
onset_control_table = arrange_tables(onset_control_feature_table, common_components, 'start')
offset_control_table = arrange_tables(offset_control_feature_table, common_components, 'end')

# add overall component score to features
if len(seizure_score_table) > 0:
    onset_seizure_table = pd.concat([onset_seizure_table, seizure_score_table], axis=1, join='inner')
    offset_seizure_table = pd.concat([offset_seizure_table, seizure_score_table], axis=1, join='inner')
    onset_control_table = pd.concat([onset_control_table, control_score_table], axis=1, join='inner')
    offset_control_table = pd.concat([offset_control_table, control_score_table], axis=1, join='inner')

print("correlation between brain scores")
if len(seizure_score_table) == 0:
    print(ss.pearsonr(onset_seizure_table['brain'], onset_seizure_table['mara']))
    print(ss.pearsonr(onset_control_table['brain'], onset_control_table['mara']))
else:
    print(ss.pearsonr(seizure_score_table['brain_both'], seizure_score_table['mara_both']))
    print(ss.pearsonr(control_score_table['brain_both'], control_score_table['mara_both']))

onset_control_table['bif_start'] = 'Control'
offset_control_table['bif_end'] = 'Control'

onset_control_table['clear_start'] = 'Control'
offset_control_table['clear_end'] = 'Control'

# define one clear
clear_ind = (offset_seizure_table['clear_end'] == 'Clear') | (onset_seizure_table['clear_start'] == 'Clear')
onset_seizure_table.loc[:, 'one_clear'] = 'Unclear'
onset_seizure_table.loc[clear_ind, 'one_clear'] = 'Clear'
offset_seizure_table.loc[:, 'one_clear'] = 'Unclear'
offset_seizure_table.loc[clear_ind, 'one_clear'] = 'Clear'
onset_control_table['one_clear'] = 'Control'
offset_control_table['one_clear'] = 'Control'
DATA_DIR = Path(DATA_DIR) / 'feature_analysis_results'
# save joined table
os.makedirs(DATA_DIR, exist_ok=True)

onset_seizure_table.reset_index(drop=False).to_csv(join(DATA_DIR, 'onset_table.csv'), index_label='index')
offset_seizure_table.reset_index(drop=False).to_csv(join(DATA_DIR, 'offset_table.csv'), index_label='index')

onset_control_table.reset_index(drop=False).to_csv(join(DATA_DIR, 'onset_control_table.csv'),
                                                   index_label='index')
offset_control_table.reset_index(drop=False).to_csv(join(DATA_DIR, 'offset_control_table.csv'),
                                                    index_label='index')


# change point analysis
def create_changepoint_summary(table, control_table, start_end):
    grouped_control = control_table.groupby([f'clear_{start_end}', 'Change point']).count()['bif_time']
    grouped_control = grouped_control / len(control_table)

    grouped = table.groupby([f'clear_{start_end}', 'Change point']).count()['bif_time']
    grouped = grouped / table.groupby([f'clear_{start_end}']).count()['bif_time']

    plot_data = pd.concat([grouped, grouped_control]) * 100
    plot_data = plot_data.reset_index(drop=False).rename(columns={'bif_time': 'percent'})
    return plot_data


def create_proportion_plot_and_save(plot_data_onset, plot_data_offset, x, y, what, subset, palette):
    fig, ax = plt.subplots(1, 2, figsize=(15, 7), sharey='row')
    sns.barplot(x=f'{x}_start', y=y,
                data=plot_data_onset[plot_data_onset[what] == subset],
                palette=palette, ax=ax[0])

    ax[0].set_title(f'Onset')
    ax[0].set_ylabel('% Detected')
    ax[0].set_xlabel('')
    ax[0].set_ylim([0, 100])

    sns.barplot(x=f'{x}_end', y=y,
                data=plot_data_offset[plot_data_offset[what] == subset],
                palette=palette, ax=ax[1])

    ax[1].set_title(f'Offset')
    ax[1].set_xlabel('')
    ax[1].set_ylabel('')
    ax[1].set_ylim([0, 100])

    labels = list(group_palette.keys())
    handles = [plt.Rectangle((0, 0), 1, 1, color=palette[label]) for label in labels]
    ax[1].legend(handles, labels, bbox_to_anchor=(1.35, 1))

    fig.suptitle('Detected change point by variance', fontsize=40)
    fig.tight_layout(pad=2)
    fig_name = f'identified_{what}_{DATA_ALIAS}'
    print(fig_name)
    os.makedirs(FIG_DIR, exist_ok=True)
    fig.savefig(os.path.join(FIG_DIR, f'{fig_name}.svg'), transparent=True, bbox_inches='tight', pad_inches=1)


plot_data_onset = create_changepoint_summary(onset_seizure_table, onset_control_table, 'start')
plot_data_offset = create_changepoint_summary(offset_seizure_table, offset_control_table, 'end')

# save the change point proportion table
pd.concat([plot_data_onset, plot_data_offset], axis=1).to_csv(
    os.path.join(DATA_DIR, 'change_point_summary_table.csv'), index_label='index')

# plot changepoint proportions
create_proportion_plot_and_save(plot_data_onset, plot_data_offset,
                                x='clear', y='percent', what='Change point', subset='Detected',
                                palette=group_palette)

# aggregate onset and offset data with control
onset_table_all = aggregate_control(onset_seizure_table, onset_control_table)
offset_table_all = aggregate_control(offset_seizure_table, offset_control_table)
drop_col = ['file_name', 'bif_start', 'bif_end', 'clear_end', 'clear_start', 'one_clear', 'Change point',
            'index.1', 'scm', 'notes', 'count', 'use_sim', 'file_name']  # 'hight_autocorr'
# score analysis - does this
score_type = 'score'
score_group_statistics = pd.DataFrame()
score_group_statistics['all control components'] = onset_control_table.describe()[score_type]
score_group_statistics['one clear'] = onset_seizure_table.loc[clear_ind, :].describe()[score_type]
score_group_statistics['none clear'] = onset_seizure_table.loc[~clear_ind, :].describe()[score_type]
score_group_statistics['all seizure components'] = onset_seizure_table.describe()[score_type]
score_group_statistics.to_csv(join(DATA_DIR, f'{score_type}_statistics.csv'), index_label='measure')


def reformat_post_hoc(post_hoc):
    for i in range(0, len(post_hoc)):
        data = post_hoc.iloc[i + 1:, i]
        data.rename(lambda x: f'{x}-{data.name}', inplace=True)
        if i == 0:
            post_hoc_pairs = data.copy()
        else:
            post_hoc_pairs = pd.concat([post_hoc_pairs, data])
    return post_hoc_pairs.to_dict()


def get_discriminative_features(table, target, drop_col, p):
    # remove rows with more than 10 0 values
    table_temp = table.loc[(table == 0).sum(axis=1) <= 10, :]
    table_temp = table_temp.replace([np.inf, -np.inf], np.nan)
    table_temp = table_temp.dropna(how="any")
    anova = pd.DataFrame(data=f_classif(table_temp.drop(columns=drop_col), table_temp[target]),
                         index=['F', 'P'],
                         columns=table_temp.drop(columns=drop_col).columns.values)

    anova = anova.T.sort_values(by="P")
    features = sorted(list(anova[anova['P'] < p].T.columns))
    return anova, features


def anova_with_post_hoc(data: pd.DataFrame, feature: str, target: str):
    lm = sfa.ols(f'{feature} ~ C({target})', data=data).fit()
    anova = sa.stats.anova_lm(lm)
    # print(anova)
    post_hoc = sp.posthoc_ttest(data, val_col=feature, group_col=target, p_adjust='holm')
    # print(post_hoc)
    return {f'{feature}-{target}_anova': {'anova': anova.to_dict(), 'post_hoc_anova': reformat_post_hoc(post_hoc)}}


def kruskal_with_post_hoc(data: pd.DataFrame, feature: str, target: str):
    data_list = [data.loc[ids, feature].values for ids in data.groupby(target).groups.values()]
    H, p = ss.kruskal(*data_list)
    # print(f'H: {H}, p: {p}')
    post_hoc = sp.posthoc_conover(data, val_col=feature, group_col=target, p_adjust='holm')
    # print(post_hoc)
    return {f'{feature}-{target}_kruskal': {'H_kruskal': H, 'p_kruskal': p,
                                            'post_hoc_kruskal': reformat_post_hoc(post_hoc)}}


# create seizure group names
all_columns = onset_table_all.columns
score_features = ['max_score', 'score', 'mara', 'brain', 'eye', 'muscle', 'chan', 'line', 'heart', 'other']

only_features = list(all_columns.drop(drop_col).values)
isi_features = list(filter(lambda x: 'isi_' in x, all_columns))
rms_features = list(filter(lambda x: 'rms_' in x, all_columns))
peak_features = list(filter(lambda x: 'peak_' in x, all_columns))
statistical_features = list(filter(lambda x: x in ['mean', 'median', 'stds', 'skewness'], all_columns))
spectral_features = list(filter(lambda x: 'spectral' in x, all_columns))
autocorr_features = list(filter(lambda x: 'autocor' in x, all_columns))
# autocorr_features.remove('height_autocorr')
entropy_features = list(filter(lambda x: 'entropy' in x, all_columns))

# group_features = list(set(score_features + spectral_features + autocorr_features + entropy_features + ['comp_num']))
# bif_features = list(
#     set(isi_features + rms_features + peak_features + autocorr_features + entropy_features + statistical_features))

group_features = list(only_features)
bif_features = list(only_features)

p_th = 0.05 / len(only_features)  # correct for the number of features tested

anova_onset_score, use_to_separate_onset = get_discriminative_features(
    onset_table_all[group_features + ['clear_start']],
    'clear_start', ['clear_start'], p_th)

anova_offset_score, use_to_separate_offset = get_discriminative_features(
    offset_table_all[group_features + ['clear_end']],
    'clear_end', ['clear_end'], p_th)

anova_onset_score_only_seizure, use_to_separate_onset_only_seizure = get_discriminative_features(
    onset_seizure_table.reset_index()[list(set(group_features + bif_features)) + ['clear_start']],
    'clear_start', ['clear_start'], p_th)

anova_offset_score_only_seizure, use_to_separate_offset_only_seizure = get_discriminative_features(
    offset_seizure_table.reset_index()[list(set(group_features + bif_features)) + ['clear_end']],
    'clear_end', ['clear_end'], p_th)

anova_onset_score_bif, use_to_separate_onset_bif = get_discriminative_features(
    onset_seizure_table.reset_index()[bif_features + ['bif_start']],
    'bif_start', ['bif_start'], p_th)

anova_offset_score_bif, use_to_separate_offset_bif = get_discriminative_features(
    offset_seizure_table.reset_index()[bif_features + ['bif_end']],
    'bif_end', ['bif_end'], p_th)

anova_onset_score.to_csv(join(DATA_DIR, 'anova_onset_group.csv'), index_label='feature')
anova_offset_score.to_csv(join(DATA_DIR, 'anova_offset_group.csv'), index_label='feature')

anova_onset_score_only_seizure.to_csv(join(DATA_DIR, 'anova_onset_group_only_seizure.csv'), index_label='feature')
anova_offset_score_only_seizure.to_csv(join(DATA_DIR, 'anova_offset_group_only_seizure.csv'), index_label='feature')

anova_onset_score_bif.to_csv(join(DATA_DIR, 'anova_onset_bif.csv'), index_label='feature')
anova_offset_score_bif.to_csv(join(DATA_DIR, 'anova_offset_bif.csv'), index_label='feature')
group_features = list(set(use_to_separate_onset) | set(use_to_separate_onset))
bif_features = list(set(use_to_separate_onset_bif) | set(use_to_separate_onset_bif))

# perform post hook tests for significant features
statistical_results_groups = dict(zip(group_features, [None] * len(group_features)))
statistical_results_between_bif_with_unclear = dict(zip(bif_features, [None] * len(bif_features)))
statistical_results_between_bif = dict(zip(bif_features, [None] * len(bif_features)))

onset_table_clear = onset_seizure_table[onset_seizure_table['clear_start'] == 'Clear']
offset_table_clear = offset_seizure_table[offset_seizure_table['clear_end'] == 'Clear']
# use two stat tests and a different package to varify and get post hoc significance
for feature in group_features:
    statistical_results_groups[feature] = []
    statistical_results_groups[feature].append(kruskal_with_post_hoc(onset_table_all, feature, 'clear_start'))
    statistical_results_groups[feature].append(kruskal_with_post_hoc(offset_table_all, feature, 'clear_end'))
    statistical_results_groups[feature].append(anova_with_post_hoc(onset_table_all, feature, 'clear_start'))
    statistical_results_groups[feature].append(anova_with_post_hoc(offset_table_all, feature, 'clear_end'))

for feature in bif_features:
    # including unclear
    statistical_results_between_bif_with_unclear[feature] = []

    statistical_results_between_bif_with_unclear[feature].append(
        kruskal_with_post_hoc(onset_seizure_table.reset_index(), feature, 'bif_start'))
    statistical_results_between_bif_with_unclear[feature].append(
        kruskal_with_post_hoc(offset_seizure_table.reset_index(), feature, 'bif_end'))
    statistical_results_between_bif_with_unclear[feature].append(
        anova_with_post_hoc(onset_seizure_table.reset_index(), feature, 'bif_start'))
    statistical_results_between_bif_with_unclear[feature].append(
        anova_with_post_hoc(offset_seizure_table.reset_index(), feature, 'bif_end'))

    # only clear
    statistical_results_between_bif[feature] = []
    statistical_results_between_bif[feature].append(
        kruskal_with_post_hoc(onset_table_clear.reset_index(), feature, 'bif_start'))
    statistical_results_between_bif[feature].append(
        kruskal_with_post_hoc(offset_table_clear.reset_index(), feature, 'bif_end'))
    statistical_results_between_bif[feature].append(
        anova_with_post_hoc(onset_table_clear.reset_index(), feature, 'bif_start'))
    statistical_results_between_bif[feature].append(
        anova_with_post_hoc(offset_table_clear.reset_index(), feature, 'bif_end'))

with open(join(DATA_DIR, 'statistical_results_groups.json'), 'w') as f:
    json.dump(statistical_results_groups, f, indent=4, sort_keys=True)

with open(join(DATA_DIR, 'statistical_results_between_bif_with_unclear.json'), 'w') as f:
    json.dump(statistical_results_between_bif_with_unclear, f, indent=4, sort_keys=True)

with open(join(DATA_DIR, 'statistical_results_between_bif.json'), 'w') as f:
    json.dump(statistical_results_between_bif, f, indent=4, sort_keys=True)

for feature in bif_features:
    plot_feature_and_save(feature, onset_seizure_table.reset_index(), offset_seizure_table.reset_index(),
                          onset_palette_only_bif, offset_palette_only_bif,
                          join(FIG_DIR, 'without_unclear'), ylim=None, error_bars=False)

for feature in bif_features:
    plot_feature_and_save(feature, onset_seizure_table.reset_index(), offset_seizure_table.reset_index(),
                          onset_palette_seizure, offset_palette_seizure,
                          join(FIG_DIR, 'with_unclear'), ylim=None, error_bars=False)

for feature in group_features:
    plot_feature_and_save(feature, onset_table_all.reset_index(), offset_table_all.reset_index(),
                          group_palette, group_palette,
                          join(FIG_DIR, 'group'), ylim=None, bif_clear='clear', error_bars=False)
