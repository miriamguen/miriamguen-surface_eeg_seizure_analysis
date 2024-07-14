""" Define the functions used in labeling_validation_analysis """

from pathlib import Path
import os
import warnings
import pickle

from tqdm import tqdm
from joblib import Parallel, delayed

import pandas as pd
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pylab

from scipy import stats
from statsmodels.stats.inter_rater import cohens_kappa
from sklearn.metrics import confusion_matrix

warnings.filterwarnings("ignore", category=DeprecationWarning)
plt.style.reload_library()
plt.style.use(["science"])
params = {
    "font.family": "sans-serif",
    "font.sans-serif": "calibri",
    "text.usetex": False,
    "axes.titlesize": 40,
    "axes.labelsize": 28,
    "axes.grid": False,
    "axes.grid.axis": "y",
    "axes.grid.which": "major",
    "axes.edgecolor": "black",
    "axes.facecolor": "white",
    "axes.linewidth": 0.4,
    "ytick.labelsize": 22,
    "xtick.labelsize": 22,
    "legend.facecolor": "white",
    "legend.fontsize": 22,
    "figure.facecolor": "white",
    "savefig.transparent": True,
    "backend": "svg",  # 'Qt5Agg' #'svg'
}
pylab.rcParams.update(params)
# N_PERMUTATIONS = 10000
# BASE_PATH = Path(os.getcwd())
# DATA_PATH = BASE_PATH / "example_data" / "labels" / "our_labels"
# RESULT_PATH = BASE_PATH / "example_data" / "labels" / "analysis_results"
# FIG_PATH = BASE_PATH / "figures" / "label analysis figures"
# FEATURE_ALIAS = "example_data"
# RESULT_PATH.mkdir(parents=True, exist_ok=True)
# FIG_PATH.mkdir(parents=True, exist_ok=True)


def permute_seizure(seizure_labels: pd.DataFrame, column: str) -> pd.Series:
    """
    Permutes the specified column of seizure labels.

    Parameters:
    seizure_labels (pd.DataFrame): DataFrame containing seizure labels.
    column (str): Column to be permuted.

    Returns:
    pd.Series: Permuted column sorted by index.
    """
    seizure_labels_ = seizure_labels.copy(deep=True)
    seizure_labels_[column] = np.random.permutation(seizure_labels_[column].values)
    return seizure_labels_.set_index(["file_name", "comp_num"])[column].sort_index()


def compute_kappa_within(
    grouped_labels: pd.core.groupby.DataFrameGroupBy, column_1: str, column_2: str
) -> float:
    """
    Computes Cohen's kappa for the given columns.

    Parameters:
    grouped_labels (pd.groupby): DataFrame containing the labels.
    column_1 (str): First column for comparison.
    column_2 (str): Second column for comparison.

    Returns:
    float: Cohen's kappa value.
    """
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    shuffled_1 = grouped_labels.apply(permute_seizure, column_1)
    shuffled_2 = grouped_labels.apply(permute_seizure, column_2)

    cm = confusion_matrix(shuffled_1, shuffled_2)
    cm_df = pd.DataFrame(
        cm, index=[i for i in set(shuffled_1)], columns=[i for i in set(shuffled_2)]
    )
    return cohens_kappa(cm_df)["kappa"]


def compute_kappa(labels: pd.DataFrame, column_1: str, column_2: str) -> float:
    """
    Computes Cohen's kappa for the given columns.

    Parameters:
    labels (pd.DataFrame): DataFrame containing the labels.
    column_1 (str): First column for comparison.
    column_2 (str): Second column for comparison.

    Returns:
    float: Cohen's kappa value.
    """
    shuffled_1 = np.random.permutation(labels[column_1].values)
    shuffled_2 = np.random.permutation(labels[column_2].values)
    cm = confusion_matrix(shuffled_1, shuffled_2)
    cm_df = pd.DataFrame(
        cm, index=[i for i in set(shuffled_1)], columns=[i for i in set(shuffled_2)]
    )
    return cohens_kappa(cm_df)["kappa"]


def simulate_kappa_within_seizure_permutation(
    labels: pd.DataFrame, column_1: str, column_2: str, grouping: str = "file_name"
) -> tuple:
    """
    Simulates Cohen's kappa by permuting labels within each seizure.

    Parameters:
    labels (pd.DataFrame): DataFrame containing the labels.
    column_1 (str): First column for comparison.
    column_2 (str): Second column for comparison.
    grouping (str): Column name used for grouping (default is 'file_name').

    Returns:
    tuple: Vector of kappa values, mean, standard error of the mean, lower and upper confidence intervals.
    """
    grouped_labels = labels.groupby(grouping, group_keys=False)
    num_processors = os.cpu_count()
    kapa_vec = Parallel(n_jobs=int(num_processors / 1.25))(
        delayed(compute_kappa_within)(grouped_labels, column_1, column_2)
        for _ in tqdm(range(N_PERMUTATIONS))
    )

    sem = stats.sem(kapa_vec)
    mean = np.mean(kapa_vec)
    confidence = 0.95  # 95% confidence interval
    n = len(kapa_vec)
    h = sem * stats.t.ppf((1 + confidence) / 2, n - 1)

    ci_lower = mean - h
    ci_upper = mean + h

    print(f"Mean: {mean}")
    print(f"95% Confidence interval: ({ci_lower}, {ci_upper})")
    return (kapa_vec, mean, sem, ci_lower, ci_upper)


def simulate_kappa_overall_seizure_permutation(
    labels: pd.DataFrame, column_1: str, column_2: str
) -> tuple:
    """
    Simulates Cohen's kappa by permuting labels overall.

    Parameters:
    labels (pd.DataFrame): DataFrame containing the labels.
    column_1 (str): First column for comparison.
    column_2 (str): Second column for comparison.

    Returns:
    tuple: Vector of kappa values, mean, standard error of the mean, lower and upper confidence intervals.
    """
    kapa_vec = Parallel(n_jobs=-1)(
        delayed(compute_kappa)(labels, column_1, column_2)
        for _ in tqdm(range(N_PERMUTATIONS))
    )

    sem = stats.sem(kapa_vec)
    mean = np.mean(kapa_vec)
    confidence = 0.95  # 95% confidence interval
    n = len(kapa_vec)
    h = sem * stats.t.ppf((1 + confidence) / 2, n - 1)

    ci_lower = mean - h
    ci_upper = mean + h

    print(f"Mean: {mean}")
    print(f"95% Confidence interval: ({ci_lower}, {ci_upper})")
    return (kapa_vec, mean, sem, ci_lower, ci_upper)


def analyze_agreement(
    data: pd.DataFrame,
    column_1: str,
    column_2: str,
    result_path: Path,
    grouping: str = "file_name",
) -> tuple:
    """
    Analyzes agreement between two sets of labels.

    Parameters:
    data (pd.DataFrame): DataFrame containing the labels.
    column_1 (str): First column for comparison.
    column_2 (str): Second column for comparison.
    result_path (Path): Path to save the result.
    grouping (str): Column name used for grouping (default is 'file_name').

    Returns:
    tuple: Detailed result dictionary and summary dictionary.
    """
    labels = list(set(data[column_1]).union(set(data[column_1])))
    cm = confusion_matrix(data[column_1], data[column_2], labels=labels)
    cm_df = pd.DataFrame(cm, index=labels, columns=labels)
    print(cm_df)

    # Calculate empirical kappa
    res_start = cohens_kappa(cm_df)
    print(res_start)

    # Perform kappa simulation permuting labels overall
    print("Performing overall permutation analysis")
    (kapa_vec_o, mean_o, sem_o, ci_lower_o, ci_upper_o) = (
        simulate_kappa_overall_seizure_permutation(
            labels=data, column_1=column_1, column_2=column_2
        )
    )
    p_value_o = np.sum(kapa_vec_o >= res_start["kappa"]) / len(kapa_vec_o)

    # Perform kappa simulation permuting labels within a seizure only
    print("Performing within seizure permutation analysis")
    (kapa_vec_w, mean_w, sem_w, ci_lower_w, ci_upper_w) = (
        simulate_kappa_within_seizure_permutation(
            labels=data, column_1=column_1, column_2=column_2, grouping=grouping
        )
    )
    p_value_w = np.sum(kapa_vec_w >= res_start["kappa"]) / len(kapa_vec_w)

    result = {
        "empirical_kappa": res_start,
        "empirical_confusion": cm_df,
        "kapa_vec_within_seizure": kapa_vec_w,
        "mean_within_seizure": mean_w,
        "p_value_within_seizure": p_value_w,
        "sem_within_seizure": sem_w,
        "ci_lower_within_seizure": ci_lower_w,
        "ci_upper_within_seizure": ci_upper_w,
        "kapa_vec_overall": kapa_vec_o,
        "mean_overall": mean_o,
        "p_value_overall": p_value_o,
        "sem_overall": sem_o,
        "ci_lower_overall": ci_lower_o,
        "ci_upper_overall": ci_upper_o,
    }

    with open(result_path, "wb") as file:
        pickle.dump(result, file)

    summary = {
        "Empirical kappa": res_start["kappa"],
        "Mean kappa within seizure permutation": mean_w,
        "P-value kappa within seizure permutation": p_value_w,
        "Mean kappa overall permutation": mean_o,
        "P-value kappa overall permutation": p_value_o,
    }

    return result, summary


def load_features(file_dir):
    # load the feature data
    file_list = os.listdir(file_dir)
    source_files = list(
        set(
            map(
                lambda x: x.split(".")[0],
                filter(lambda x: x.endswith(".csv"), file_list),
            )
        )
    )

    onset_features = []
    offset_features = []

    for file in source_files:
        onset = pd.read_csv(
            file_dir / f"{file}.mat_start_bif_variance_both.csv"
        ).set_index(["file_name", "comp_num"])

        offset = pd.read_csv(
            file_dir / f"{file}.mat_end_bif_variance_both.csv"
        ).set_index(["file_name", "comp_num"])

        general = pd.read_csv(
            file_dir / f"{file}.mat_general_score_variance_both.csv"
        ).set_index(["file_name", "comp_num"])
        general = general.rename(columns=lambda x: f"{x}_general")

        onset_features.append(pd.concat([onset, general], axis=1))
        offset_features.append(pd.concat([offset, general], axis=1))

    onset_features_df = pd.concat(onset_features)
    offset_features_df = pd.concat(offset_features)

    return onset_features_df, offset_features_df


def return_seizure_dynamotype(df_input, start_end):
    # in seizures with multiple components select the component with either both bif or the highest score for analysis
    df = df_input[df_input[f"clear_bifurcation_{start_end}"]]
    df["overall_brain"] = df[["brain_general", "brain", "mara", "mara_general"]].max(
        axis=1
    )
    number_of_components = len(df)
    if number_of_components == 1:
        df["number_of_components"] = 1
        return df
    else:
        # prefer components with clear changepoint and max score
        clear_time_ind = df[df["bif_time"] > 0].index.values
        good_brain_ind = df[df["overall_brain"] > 0.6].index.values

        both = list(set(good_brain_ind) & set(clear_time_ind))

        # if there are components with good score & clear onset select the edge time
        if len(both) > 0:
            if start_end == "start":
                ind = df.loc[both, "bif_time"].idxmin()
            else:
                ind = df.loc[both, "bif_time"].idxmax()

        # if the components have no objective onset time select the highest brain score
        elif len(good_brain_ind) > 0:
            ind = df.loc[good_brain_ind, "overall_brain"].idxmax()

        # if all the components have a relatively bad brain score but do have clear time
        elif len(clear_time_ind) > 0:
            if start_end == "start":
                ind = df.loc[clear_time_ind, "bif_time"].idxmin()
            else:
                ind = df.loc[clear_time_ind, "bif_time"].idxmax()
        # if no clear onset select the brain
        else:
            ind = df["brain_general"].idxmax()

        return_df = df.loc[[ind], :]
        return_df["number_of_components"] = number_of_components
    return return_df


def plot_the_overall_distribution_pie(
    labels_start: pd.DataFrame,
    labels_end: pd.DataFrame,
    bif_start_colors: dict,
    bif_end_colors: dict,
    fig_path: Path,
) -> tuple:
    """
    Plots pie charts for the overall distribution of bifurcation types at onset and offset.

    Parameters:
    labels start (pd.DataFrame): DataFrame containing bifurcation labels.
    labels end (pd.DataFrame): DataFrame containing bifurcation labels.
    bif_start_colors (dict): Dictionary of colors for bifurcation start types.
    bif_end_colors (dict): Dictionary of colors for bifurcation end types.

    Returns:
    tuple: A tuple containing the figure, bifurcation start counts, and bifurcation end counts.
    """
    # Filter out 'Not Clear' bifurcation labels

    # Count bifurcation types
    bif_start_count = (
        labels_start["bifurcation_start"]
        .value_counts()
        .reindex(bif_start_colors.keys(), fill_value=0)
    )
    bif_end_count = (
        labels_end["bifurcation_end"]
        .value_counts()
        .reindex(bif_end_colors.keys(), fill_value=0)
    )

    # Set up the figure and title
    fig, axs = plt.subplots(1, 2, figsize=(15, 8), dpi=300)
    # fig.suptitle(title, fontsize=40, color="navy")
    textprops = {"fontsize": 30}

    # Plot start bifurcations
    axs[0].pie(
        bif_start_count,
        colors=[bif_start_colors.get(x) for x in bif_start_count.index],
        autopct="%1.1f%%",
        pctdistance=1.2,
        textprops=textprops,
        wedgeprops={"linewidth": 2, "edgecolor": "black", "alpha": 0.5},
        startangle=109.4,
        radius=1,
        counterclock=True,
    )
    axs[0].legend(
        labels=bif_start_count.index, bbox_to_anchor=(-0.3, 0.7, 0.5, 0.5), fontsize=30
    )
    axs[0].set_title(f"Onset (n ={len(labels_start)})", fontsize=36, y=1.20, x=0.5)

    # Plot end bifurcations
    axs[1].pie(
        bif_end_count,
        colors=[bif_end_colors.get(x) for x in bif_end_count.index],
        autopct="%1.1f%%",
        pctdistance=1.2,
        textprops=textprops,
        wedgeprops={"linewidth": 1.5, "edgecolor": "black", "alpha": 0.5},
        startangle=38.7,
        radius=1,
        counterclock=False,
    )
    axs[1].legend(
        labels=bif_end_count.index, fontsize=30, bbox_to_anchor=(1, 0.7, 0.5, 0.5)
    )
    axs[1].set_title(f"Offset (n={len(labels_end)})", fontsize=36, y=1.20, x=0.5)
    # fig.tight_layout(pad=1)
    # Save the figure
    fig.savefig(
        fig_path / "bifurcation_distribution.svg",
        transparent=True,
        bbox_inches="tight",
        pad_inches=1,
    )

    return fig, bif_start_count, bif_end_count


def get_proportions(df_orig, column1, column2, save_path):
    # Ensure columns are categorical
    df = df_orig.copy(deep=True)
    df = df[~(df[column1] == "unclear")]
    df = df[~(df[column2] == "unclear")]
    df[column1] = df[column1].astype("category")
    df[column2] = df[column2].astype("category")

    # Calculate the occurrence proportions
    proportions = (
        df.groupby([column1, column2], group_keys=False, observed=False)
        .size()
        .groupby(level=0, group_keys=False, observed=False)
        .apply(lambda x: x / float(x.sum()))
        .reset_index(name="proportion")
    )
    pivot_table = proportions.pivot(index=column1, columns=column2, values="proportion")
    # Plotting the heatmap
    fig = plt.figure(figsize=(6, 5))
    sns.heatmap(
        pivot_table,
        vmin=0,
        vmax=1,
        annot=True,
        cmap="cividis_r",
        cbar_kws={"label": "Proportion"},
        annot_kws={"size": 16},
        linewidths=0,
        linecolor=None,
        alpha=0.5,
    )
    plt.xlabel(column2.replace("_", " ").capitalize())
    plt.ylabel(column1.replace("_", " ").capitalize())
    # plt.title(f'Proportions of {column2.capitalize()} by {column1.capitalize()}')
    plt.grid(False)
    # Save the figure
    plt.savefig(save_path, format="svg", dpi=300)
    plt.close()

    # Pivot the data for plotting
    return pivot_table, fig
