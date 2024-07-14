""" This scrip will run the labeling analysis, inter-rater agreement and figure generation presented in the ms
    The code relays on the functions defined in analysis_utils.py
"""

from pathlib import Path
import os
import warnings

import pandas as pd
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns


from analysis_utils import *

warnings.filterwarnings("ignore", category=DeprecationWarning)


# Arrange the data
BASE_PATH = Path(os.getcwd())
DATA_PATH = BASE_PATH / "example_data" / "labels" / "our_labels"
RESULT_PATH = BASE_PATH / "example_data" / "labels" / "analysis_results"
FIG_PATH = BASE_PATH / "figures" / "label_analysis_figures"
FEATURE_ALIAS = "example_data"
RESULT_PATH.mkdir(parents=True, exist_ok=True)
FIG_PATH.mkdir(parents=True, exist_ok=True)

N_PERMUTATIONS = 10000
results = {
    "reviewer 1": {},
    "reviewer 2": {},
    "shared clear transition": {},
    "shared labeled bifurcation": {},
    "seizures with multiple bifurcations": {},
}


# %% get all the files and arrange the data, store general measures
labels_1 = pd.read_csv(DATA_PATH / "reviewer_1_labels.csv", index_col="index")
labels_2 = pd.read_csv(DATA_PATH / "reviewer_2_labels.csv", index_col="index")

file_list = []
component = []
for file in set(labels_1.reset_index(drop=False).file_name):
    for i in range(1, 9):
        file_list.append(file)
        component.append(i)

unclear_seizures_1 = labels_1[labels_1.comp_num == 0]
unclear_seizures_2 = labels_2[labels_2.comp_num == 0]

labels_1 = labels_1[labels_1.comp_num > 0].set_index(["file_name", "comp_num"])
labels_2 = labels_2[labels_2.comp_num > 0].set_index(["file_name", "comp_num"])

labels_1_ = labels_1.rename(columns=lambda x: f"{x}_1")
labels_2_ = labels_2.rename(columns=lambda x: f"{x}_2")

clear_1_start = labels_1_[~(labels_1_.loc[:, "bif_start_1"] == "not_clear")]
clear_1_end = labels_1_[~(labels_1_.loc[:, "bif_end_1"] == "not_clear")]
print("Reviewer 1 selected:\n")
print(f" Onset: {len(clear_1_start)}")
print(f" Offset: {len(clear_1_end)}\n")

results["reviewer 1"]["components with onset"] = len(clear_1_start)
results["reviewer 1"]["components with offset"] = len(clear_1_end)

clear_2_start = labels_2_[~(labels_2_.loc[:, "bif_start_2"] == "not_clear")]
clear_2_end = labels_2_[~(labels_2_.loc[:, "bif_end_2"] == "not_clear")]
print("Reviewer 2 selected:\n")
print(f" Onset: {len(clear_2_start)}")
print(f" Offset: {len(clear_2_end)}\n")
results["reviewer 2"]["components with onset"] = len(clear_2_start)
results["reviewer 2"]["components with offset"] = len(clear_2_end)

shared_onset = pd.concat([clear_1_start, clear_2_start], join="inner", axis=1)
results["shared clear transition"]["components with onset"] = len(shared_onset)
results["shared clear transition"]["seizures with onset"] = len(
    set(shared_onset.reset_index()["file_name"])
)
results["shared clear transition"]["mean number of seizure components onset"] = (
    shared_onset.reset_index().groupby(["file_name"]).count().clear_start_1.mean()
)
shared_onset_bif = shared_onset[
    shared_onset["bif_start_1"] == shared_onset["bif_start_2"]
]
results["shared labeled bifurcation"]["components with onset"] = len(shared_onset_bif)
shared_onset_bif_seizure = shared_onset_bif.reset_index().groupby(["file_name"]).sum()
results["shared labeled bifurcation"]["seizures with onset"] = len(
    shared_onset_bif_seizure
)
results["shared labeled bifurcation"][
    "mean number of seizure components onset"
] = shared_onset_bif_seizure.clear_start_1.mean()

print("\nCommon components with clear onset transition", len(shared_onset))
print(
    "Common seizures with clear onset transition",
    len(set(shared_onset.reset_index()["file_name"])),
)
print("\nSame onset bifurcation", len(shared_onset_bif))
print(
    "Seizures with same onset bifurcation",
    len(set(shared_onset_bif.reset_index()["file_name"])),
)
shared_onset_bif.to_csv(RESULT_PATH / "agreed_onset.csv", index_label="index")

shared_offset = pd.concat([clear_1_end, clear_2_end], join="inner", axis=1)

results["shared clear transition"]["components with offset"] = len(shared_offset)
results["shared clear transition"]["seizures with offset"] = len(
    set(shared_offset.reset_index()["file_name"])
)
results["shared clear transition"]["mean number of seizure components offset"] = (
    shared_offset.reset_index().groupby(["file_name"]).count().clear_end_1.mean()
)
shared_offset_bif = shared_offset[
    shared_offset["bif_end_1"] == shared_offset["bif_end_2"]
]
results["shared labeled bifurcation"]["components with offset"] = len(shared_offset_bif)
shared_offset_bif_seizure = shared_offset_bif.reset_index().groupby(["file_name"]).sum()
results["shared labeled bifurcation"]["seizures with offset"] = len(
    shared_offset_bif_seizure
)
results["shared labeled bifurcation"][
    "mean number of seizure components offset"
] = shared_offset_bif_seizure.clear_end_1.mean()

print("\nCommon components with clear offset transition", len(shared_offset))
print(
    "Common seizures with clear offset transition",
    len(set(shared_offset.reset_index()["file_name"])),
)
print("\nSame offset bifurcation", len(shared_offset_bif))
print(
    "Seizures with same offset bifurcation",
    len(set(shared_offset_bif.reset_index()["file_name"])),
)
shared_offset_bif.to_csv(RESULT_PATH / "agreed_offset.csv", index_label="index")

# %% Create a dataframe with all analyzed components

all_labels = pd.DataFrame(
    columns=["label_1_start", "label_2_start", "label_1_end", "label_2_end"],
    index=pd.MultiIndex.from_arrays(
        [file_list, component], names=("file_name", "comp_num")
    ),
)

all_labels = all_labels.assign(
    label_1_start="not_clear",
    label_2_start="not_clear",
    label_1_end="not_clear",
    label_2_end="not_clear",
)

all_labels.loc[clear_1_start.index, "label_1_start"] = clear_1_start["bif_start_1"]
all_labels.loc[clear_2_start.index, "label_2_start"] = clear_2_start["bif_start_2"]
all_labels.loc[clear_1_end.index, "label_1_end"] = clear_1_end["bif_end_1"]
all_labels.loc[clear_2_end.index, "label_2_end"] = clear_2_end["bif_end_2"]

is_clear = lambda x: "not_clear" if x == "not_clear" else "clear"

all_labels["clear_1_start"] = all_labels["label_1_start"].apply(is_clear)
all_labels["clear_2_start"] = all_labels["label_2_start"].apply(is_clear)
all_labels["clear_1_end"] = all_labels["label_1_end"].apply(is_clear)
all_labels["clear_2_end"] = all_labels["label_2_end"].apply(is_clear)

all_labels = all_labels.reset_index()

results["reviewer 1"]["components reviewed overall"] = len(all_labels)
results["reviewer 2"]["components reviewed overall"] = len(all_labels)

# %% Components by seizure evaluated per reviewer
components_per_seizure_1 = all_labels.groupby("file_name", group_keys=False).apply(
    lambda x: sum(x["clear_1_start"] == "clear")
)
components_per_seizure_1 = components_per_seizure_1[components_per_seizure_1 > 0]
results["reviewer 1"]["seizures with onset"] = len(components_per_seizure_1)
results["reviewer 1"][
    "mean number of seizure components onset"
] = components_per_seizure_1.mean()

components_per_seizure_2 = all_labels.groupby("file_name", group_keys=False).apply(
    lambda x: sum(x["clear_2_start"] == "clear")
)
components_per_seizure_2 = components_per_seizure_2[components_per_seizure_2 > 0]
results["reviewer 2"]["seizures with onset"] = len(components_per_seizure_2)
results["reviewer 2"][
    "mean number of seizure components onset"
] = components_per_seizure_2.mean()

components_per_seizure_1 = all_labels.groupby("file_name", group_keys=False).apply(
    lambda x: sum(x["clear_1_end"] == "clear")
)
components_per_seizure_1 = components_per_seizure_1[components_per_seizure_1 > 0]
results["reviewer 1"]["seizures with offset"] = len(components_per_seizure_1)
results["reviewer 1"][
    "mean number of seizure components offset"
] = components_per_seizure_1.mean()

components_per_seizure_2 = all_labels.groupby("file_name", group_keys=False).apply(
    lambda x: sum(x["clear_2_end"] == "clear")
)
components_per_seizure_2 = components_per_seizure_2[components_per_seizure_2 > 0]
results["reviewer 2"]["seizures with offset"] = len(components_per_seizure_2)
results["reviewer 2"][
    "mean number of seizure components offset"
] = components_per_seizure_2.mean()

pd.DataFrame(results).sort_index().to_csv(
    RESULT_PATH / "agreement_numbers.csv", index_label="index"
)

# %% Estimate inter rater agreement

# agreement_analysis_stats = {}

# # estimate agreement on clear transition identification
# _, agreement_analysis_stats["All clear transition onset"] = analyze_agreement(
#     data=all_labels,
#     column_1="clear_1_start",
#     column_2="clear_2_start",
#     result_path=RESULT_PATH / "all_clear_transition_onset_permutation_artifacts.pkl",
#     grouping="file_name",
# )

# _, agreement_analysis_stats["All clear transition offset"] = analyze_agreement(
#     data=all_labels,
#     column_1="clear_1_end",
#     column_2="clear_2_end",
#     result_path=RESULT_PATH / "all_clear_transition_offset_permutation_artifacts.pkl",
#     grouping="file_name",
# )

# # estimate agreement on bifurcation labeling over all components
# _, agreement_analysis_stats["All labeled bifurcation onset"] = analyze_agreement(
#     data=all_labels,
#     column_1="label_1_start",
#     column_2="label_2_start",
#     result_path=RESULT_PATH / "all_bifurcation_onset_permutation_artifacts.pkl",
#     grouping="file_name",
# )

# _, agreement_analysis_stats["All labeled bifurcation offset"] = analyze_agreement(
#     data=all_labels,
#     column_1="label_1_end",
#     column_2="label_2_end",
#     result_path=RESULT_PATH / "all_bifurcation_offset_permutation_artifacts.pkl",
#     grouping="file_name",
# )

# # estimate agreement on bifurcation labeling over mutually selected clear transitions
# _, agreement_analysis_stats["Shared labeled bifurcation onset"] = analyze_agreement(
#     data=shared_onset.reset_index(),
#     column_1="bif_start_1",
#     column_2="bif_start_2",
#     result_path=RESULT_PATH
#     / "shared_transition_bifurcation_onset_permutation_artifacts.pkl",
#     grouping="file_name",
# )

# _, agreement_analysis_stats["Shared labeled bifurcation offset"] = analyze_agreement(
#     data=shared_offset.reset_index(),
#     column_1="bif_end_1",
#     column_2="bif_end_2",
#     result_path=RESULT_PATH
#     / "shared_transition_bifurcation_offset_permutation_artifacts.pkl",
#     grouping="file_name",
# )

# flattened_data = []
# for key, value in agreement_analysis_stats.items():
#     flat_dict = {"Category": key}
#     flat_dict.update(value)
#     flattened_data.append(flat_dict)

# pd.DataFrame(flattened_data).to_csv(RESULT_PATH / "agreement_statistics.csv")

# %% Load the automated features for component selection plot the bifurcation distribution plots + score + change point analysis


onset_features_df, offset_features_df = load_features(
    BASE_PATH / "example_data" / "features" / "seizure"
)


all_labels = all_labels.set_index(["file_name", "comp_num"])
all_labels["clear_start"] = np.logical_and(
    all_labels["clear_1_start"] == "clear", all_labels["clear_2_start"] == "clear"
)
all_labels["clear_bifurcation_start"] = np.logical_and(
    all_labels["clear_start"],
    all_labels["label_1_start"] == all_labels["label_2_start"],
)
all_labels["bifurcation_start"] = "not_clear"
all_labels.loc[all_labels["clear_bifurcation_start"], "bifurcation_start"] = (
    all_labels.loc[all_labels["clear_bifurcation_start"], "label_1_start"]
)


all_labels["clear_end"] = np.logical_and(
    all_labels["clear_1_end"] == "clear", all_labels["clear_2_end"] == "clear"
)
all_labels["clear_bifurcation_end"] = np.logical_and(
    all_labels["clear_end"], all_labels["label_1_end"] == all_labels["label_2_end"]
)
all_labels["bifurcation_end"] = "not_clear"
all_labels.loc[all_labels["clear_bifurcation_end"], "bifurcation_end"] = all_labels.loc[
    all_labels["clear_bifurcation_end"], "label_1_end"
]
# save for other analysis
all_labels.to_csv(RESULT_PATH / "all_labels.csv", index_label="index")


# analyze onset
onset_features_df = pd.concat(
    [
        all_labels[["clear_bifurcation_start", "clear_start", "bifurcation_start"]],
        onset_features_df,
    ],
    axis=1,
)


onset_features_df.reset_index().to_csv(
    RESULT_PATH / "all_onset_features_for_classification.csv", index_label="index"
)


# before selecting the best component record all existing
# bifurcations per seizure
all_bif_types_clear_onset = (
    onset_features_df[onset_features_df["clear_bifurcation_start"]]
    .groupby("file_name", group_keys=False)
    .apply(lambda df: set(df["bifurcation_start"].values))
)

n_bif_onset = all_bif_types_clear_onset.apply(len)
results["seizures with multiple bifurcations"]["seizures with onset"] = sum(
    n_bif_onset > 1
)


seizure_onset = (
    onset_features_df[all_labels["clear_bifurcation_start"]]
    .reset_index()
    .groupby("file_name", group_keys=False)
    .apply(return_seizure_dynamotype, "start")
    .set_index(["file_name"])
)

seizure_onset["onset_bif_types"] = all_bif_types_clear_onset
seizure_onset["onset_bif_type_number"] = n_bif_onset

onset_bif = seizure_onset["bifurcation_start"].value_counts()

# analyze offset
offset_features_df = pd.concat(
    [
        all_labels[["clear_bifurcation_end", "clear_end", "bifurcation_end"]],
        offset_features_df,
    ],
    axis=1,
)

offset_features_df.reset_index().to_csv(
    RESULT_PATH / "all_offset_features_for_classification.csv", index_label="index"
)


# before selecting the best component record all existing bifurcations per seizure
all_bif_types_clear_offset = (
    offset_features_df[offset_features_df["clear_bifurcation_end"]]
    .groupby("file_name", group_keys=False)
    .apply(lambda df: set(df["bifurcation_end"].values))
)

n_bif_offset = all_bif_types_clear_offset.apply(len)
results["seizures with multiple bifurcations"]["seizures with offset"] = sum(
    n_bif_offset > 1
)


seizure_offset = (
    offset_features_df[all_labels["clear_bifurcation_end"]]
    .reset_index()
    .groupby("file_name", group_keys=False)
    .apply(return_seizure_dynamotype, "end")
    .set_index(["file_name"])
)

seizure_offset["offset_bif_types"] = all_bif_types_clear_offset
seizure_offset["offset_bif_type_number"] = n_bif_offset

# resave the seizures with multiple bifurcations
pd.DataFrame(results).sort_index().to_csv(
    RESULT_PATH / "agreement_numbers.csv", index_label="index"
)


onset_offset = pd.concat(
    [seizure_onset["bifurcation_start"], seizure_offset["bifurcation_end"]],
    join="outer",
    axis=1,
).fillna("Not clear")

dynamotypes = pd.crosstab(
    onset_offset["bifurcation_start"], onset_offset["bifurcation_end"], margins=True
)


dynamotypes["Clear"] = dynamotypes["All"] - dynamotypes["Not clear"]
dynamotypes.loc["Clear", :] = (
    dynamotypes.loc["All", :] - dynamotypes.loc["Not clear", :]
)


dynamotypes.loc[
    ["SN/SubH", "SupH", "SNIC", "Clear", "Not clear", "All"],
    ["FLC", "SH/SNIC", "SupH", "Clear", "Not clear", "All"],
].reset_index().rename(
    columns={"bifurcation_start": "Bifurcation start / Bifurcation end"}
).to_csv(
    RESULT_PATH / "dynamotype_table.csv", index=False, header=True
)


# %% Create the bifurcation distribution pie plot


bif_start_colors = {
    "SN/SubH": "royalblue",
    "SupH": "forestgreen",
    "SNIC": "gold",
}  # , 'Multiple': 'maroon'}
bif_end_colors = {
    "FLC": "cyan",
    "SH/SNIC": "darkorange",
    "SupH": "forestgreen",
}  # , 'Multiple': 'maroon'}

onset_class_ = list(bif_start_colors.keys())
offset_class_ = list(bif_end_colors.keys())


# Example usage:
# Define your DataFrame `labels` and color dictionaries `bif_start_colors` and `bif_end_colors`
# fig, bif_start_count, bif_end_count = plot_the_overall_distribution_pie(labels, "Bifurcation Distribution", bif_start_colors, bif_end_colors)


fig, onset, offset = plot_the_overall_distribution_pie(
    seizure_onset, seizure_offset, bif_start_colors, bif_end_colors, fig_path=FIG_PATH
)

# %% Create the vigilance bifurcation confusion matrix plot

metadata = pd.read_csv(
    BASE_PATH / "example_data" / "sorce_data" / "metadata.csv", index_col=0
).set_index("file_name")


seizure_onset = pd.concat([seizure_onset, metadata], axis=1, join="inner")
# save data for clinical factors analysis
seizure_onset.to_csv(RESULT_PATH / "agreed_onset.csv", index_label="index")
seizure_offset = pd.concat([seizure_offset, metadata], axis=1, join="inner")
# save data for clinical factors analysis
seizure_offset.to_csv(RESULT_PATH / "agreed_offset.csv", index_label="index")


# %% plot the vigilance, classification confusion matrix
prop, fig = get_proportions(
    seizure_onset, "vigilance", "bifurcation_start", RESULT_PATH / "vig_prop_onset.svg"
)

fig.savefig(FIG_PATH / "vigilance_cm_onset.svg", dpi=300)

prop, fig = get_proportions(
    seizure_offset, "vigilance", "bifurcation_end", RESULT_PATH / "vig_prop_offset.svg"
)
fig.savefig(FIG_PATH / "vigilance_cm_offset.svg", dpi=300)



prop, fig = get_proportions(
    seizure_onset, "classification", "bifurcation_start", RESULT_PATH / "vig_prop_onset.svg"
)

fig.savefig(FIG_PATH / "classification_cm_onset.svg", dpi=300)

prop, fig = get_proportions(
    seizure_offset, "classification", "bifurcation_end", RESULT_PATH / "vig_prop_offset.svg"
)
fig.savefig(FIG_PATH / "classification_cm_offset.svg", dpi=300)

# %% Plot component characteristics component level
# %% Objective change point detection proportion plot
bif_start_colors = {
    "Control": "gray",
    "Unclear": "navy",
    "SN/SubH": "royalblue",
    "SupH": "forestgreen",
    "SNIC": "gold",
}  # , 'Multiple': 'maroon'}
bif_end_colors = {
    "Control": "gray",
    "Unclear": "navy",
    "FLC": "cyan",
    "SH/SNIC": "darkorange",
    "SupH": "forestgreen",
}  # , 'Multiple': 'maroon'}


control_onset_features_df, control_offset_features_df = load_features(
    BASE_PATH / "example_data" / "features" / "control"
)

onset_features_df["clear_change_point"] = (onset_features_df["bif_time"] > 0).astype(
    int
)
control_onset_features_df["clear_change_point"] = (
    control_onset_features_df["bif_time"] > 0
).astype(int)

change_point_ratio_by_bif_onset = onset_features_df.groupby(["bifurcation_start"])[
    "clear_change_point"
].mean()
change_point_ratio_by_bif_onset["Control"] = control_onset_features_df[
    "clear_change_point"
].mean()
change_point_ratio_by_bif_onset = change_point_ratio_by_bif_onset.rename(
    lambda x: x.replace("not_clear", "Unclear")
).sort_values()

offset_features_df["clear_change_point"] = (offset_features_df["bif_time"] > 0).astype(
    int
)
control_offset_features_df["clear_change_point"] = (
    control_offset_features_df["bif_time"] > 0
).astype(int)

change_point_ratio_by_bif_offset = offset_features_df.groupby(["bifurcation_end"])[
    "clear_change_point"
].mean()
change_point_ratio_by_bif_offset["Control"] = control_offset_features_df[
    "clear_change_point"
].mean()
change_point_ratio_by_bif_offset = change_point_ratio_by_bif_offset.rename(
    lambda x: x.replace("not_clear", "Unclear")
).sort_values()

# Create the plot
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(15, 6), sharey="row")

# Onset plot
sns.barplot(
    x=change_point_ratio_by_bif_onset.index,
    y=change_point_ratio_by_bif_onset.values,
    hue=change_point_ratio_by_bif_onset.index,
    palette=bif_start_colors,
    ax=axes[0],
    alpha=0.5,
)

axes[0].set_title("Onset", fontsize=28)
axes[0].set_ylabel("Detected change-point ratio")
axes[0].set_xlabel("")
# Offset plot

sns.barplot(
    x=change_point_ratio_by_bif_offset.index,
    y=change_point_ratio_by_bif_offset.values,
    hue=change_point_ratio_by_bif_offset.index,
    palette=bif_end_colors,
    ax=axes[1],
    alpha=0.5,
)

axes[1].set_title("Offset", fontsize=28)
axes[1].set_xlabel("")
plt.tight_layout()

fig.savefig(FIG_PATH / "change_point_detected_proportion.svg", dpi=300)
# %% Brain score comparing bifurcation labels and unclear, same for eye score
control_onset_features_df["bifurcation_start"] = "Control"

onset_control_seizure = pd.concat(
    [
        control_onset_features_df[["bifurcation_start", "score"]],
        onset_features_df[["bifurcation_start", "score"]],
    ],
    axis=0,
).reset_index()

onset_control_seizure["bifurcation_start"] = onset_control_seizure[
    "bifurcation_start"
].apply(lambda x: x.replace("not_clear", "Unclear"))


control_offset_features_df["bifurcation_end"] = "Control"

offset_control_seizure = pd.concat(
    [
        control_offset_features_df[["bifurcation_end", "score"]],
        offset_features_df[["bifurcation_end", "score"]],
    ],
    axis=0,
).reset_index()

offset_control_seizure["bifurcation_end"] = offset_control_seizure[
    "bifurcation_end"
].apply(lambda x: x.replace("not_clear", "Unclear"))


# Create the plot
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(15, 6), sharey="row")

# Onset plot
sns.violinplot(
    data=onset_control_seizure,
    x="bifurcation_start",
    y="score",
    hue=onset_control_seizure["bifurcation_start"],
    palette=bif_start_colors,
    ax=axes[0],
    order=["Control", "Unclear", "SNIC", "SN/SubH", "SupH"],
)
for violin in axes[0].collections:
    violin.set_alpha(0.5)

axes[0].set_title("Onset", fontsize=28)
axes[0].set_ylabel("Brain score distribution")
axes[0].set_xlabel("")
# Offset plot

sns.violinplot(
    data=offset_control_seizure,
    x="bifurcation_end",
    y="score",
    hue=offset_control_seizure["bifurcation_end"],
    palette=bif_end_colors,
    ax=axes[1],
    order=["Control", "Unclear", "SH/SNIC", "FLC", "SupH"],
)
for violin in axes[1].collections:
    violin.set_alpha(0.5)

axes[1].set_title("Offset", fontsize=28)
axes[1].set_xlabel("")
plt.tight_layout()

fig.savefig(FIG_PATH / "brain_score_by_category_proportion.svg", dpi=300)


# %% Plot detectability proportion seizure level

metadata["clear_bifurcation_start"] = False
metadata.loc[seizure_onset.index, "clear_bifurcation_start"] = True
metadata["bifurcation_start"] = None
metadata.loc[seizure_onset.index, "bifurcation_start"] = seizure_onset[
    "bifurcation_start"
]
metadata["clear_bifurcation_end"] = False
metadata.loc[seizure_offset.index, "clear_bifurcation_end"] = True
metadata["bifurcation_end"] = None
metadata.loc[seizure_offset.index, "bifurcation_end"] = seizure_offset[
    "bifurcation_end"
]


# plot detection rate by labeled onset location
import matplotlib as mpl

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(15, 6), sharey="row")

lobe_detection_ratio_onset = (
    metadata.groupby("lobe")["clear_bifurcation_start"]
    .mean()
    .sort_values()
    .rename(lambda x: x.capitalize())
)
lobe_detection_ratio_offset = (
    metadata.groupby("lobe")["clear_bifurcation_end"]
    .mean()
    .sort_values()
    .rename(lambda x: x.capitalize())
)
# Onset plot
sns.barplot(
    x=lobe_detection_ratio_onset.index,
    y=lobe_detection_ratio_onset.values,
    hue=lobe_detection_ratio_onset.values,
    palette="bone",  #'gist_gray_r',
    ax=axes[0],
    alpha=0.5,
    legend=False,
    hue_norm=mpl.colors.Normalize(vmin=0, vmax=1),
)

axes[0].set_title("Onset", fontsize=28)
axes[0].set_ylabel("Clear bifurcation ratio")
axes[0].set_xlabel("")
axes[0].tick_params(axis="x", rotation=20)
# Offset plot

sns.barplot(
    x=lobe_detection_ratio_offset.index,
    y=lobe_detection_ratio_offset.values,
    hue=lobe_detection_ratio_offset.values,
    palette="bone",  # 'gist_gray_r',
    ax=axes[1],
    alpha=0.5,
    legend=False,
    hue_norm=mpl.colors.Normalize(vmin=0, vmax=1),
)

axes[1].set_title("Offset", fontsize=28)
axes[1].set_xlabel("")
axes[1].tick_params(axis="x", rotation=20)
plt.tight_layout()

fig.savefig(FIG_PATH / "detection_rate_lobe.svg", dpi=300)


###  plot

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(15, 6), sharey="row")

lobe_detection_ratio_onset = (
    metadata.groupby("classification")["clear_bifurcation_start"]
    .mean()
    .sort_values()
    .rename({"SG": "FBTCS", "CP": "FIAS", "SP": "FAS"})
)
lobe_detection_ratio_offset = (
    metadata.groupby("classification")["clear_bifurcation_end"]
    .mean()
    .sort_values()
    .rename({"SG": "FBTCS", "CP": "FIAS", "SP": "FAS"})
)
# Onset plot for classification
sns.barplot(
    x=lobe_detection_ratio_onset.index,
    y=lobe_detection_ratio_onset.values,
    hue=lobe_detection_ratio_onset.values,
    palette="bone",  #'gist_gray_r',
    ax=axes[0],
    alpha=0.5,
    legend=False,
    hue_norm=mpl.colors.Normalize(vmin=0, vmax=1),
)

axes[0].set_title("Onset", fontsize=28)
axes[0].set_ylabel("Clear bifurcation ratio")
axes[0].set_xlabel("")
axes[0].tick_params(axis="x", rotation=0)
# Offset plot

sns.barplot(
    x=lobe_detection_ratio_offset.index,
    y=lobe_detection_ratio_offset.values,
    hue=lobe_detection_ratio_offset.values,
    palette="bone",  # 'gist_gray_r',
    ax=axes[1],
    alpha=0.5,
    legend=False,
    hue_norm=mpl.colors.Normalize(vmin=0, vmax=1),
)

axes[1].set_title("Offset", fontsize=28)
axes[1].set_xlabel("")
axes[1].tick_params(axis="x", rotation=0)
plt.tight_layout()

fig.savefig(FIG_PATH / "detection_rate_classification.svg", dpi=300)

# plot for vigilance


fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(15, 6), sharey="row")

lobe_detection_ratio_onset = (
    metadata.groupby("vigilance")["clear_bifurcation_start"].mean().sort_values()
)
lobe_detection_ratio_offset = (
    metadata.groupby("vigilance")["clear_bifurcation_end"].mean().sort_values()
)
# Onset plot
sns.barplot(
    x=lobe_detection_ratio_onset.index,
    y=lobe_detection_ratio_onset.values,
    hue=lobe_detection_ratio_onset.values,
    palette="bone",  #'gist_gray_r',
    ax=axes[0],
    alpha=0.5,
    legend=False,
    hue_norm=mpl.colors.Normalize(vmin=0, vmax=1),
)

axes[0].set_title("Onset", fontsize=28)
axes[0].set_ylabel("Clear bifurcation ratio")
axes[0].set_xlabel("")
axes[0].tick_params(axis="x", rotation=90)
# Offset plot

sns.barplot(
    x=lobe_detection_ratio_offset.index,
    y=lobe_detection_ratio_offset.values,
    hue=lobe_detection_ratio_offset.values,
    palette="bone",  # 'gist_gray_r',
    ax=axes[1],
    alpha=0.5,
    legend=False,
    hue_norm=mpl.colors.Normalize(vmin=0, vmax=1),
)

axes[1].set_title("Offset", fontsize=28)
axes[1].set_xlabel("")
axes[1].tick_params(axis="x", rotation=90)
plt.tight_layout()

fig.savefig(FIG_PATH / "detection_rate_vigilance.svg", dpi=300)


# plot the number of bifurcations by patient vs number of seizures per patient

offset_counts = metadata.groupby("patient_id")["bifurcation_end"].apply(
    lambda x: (
        x.value_counts()
        if len(x.value_counts()) > 0
        else pd.Series(name="count", index=["Unclear"], data=[None])
    )
)
offset_counts = offset_counts.unstack()
offset_counts["Offset bifurcation types"] = (offset_counts > 0).sum(axis=1)
offset_counts["Clear offset bifurcations"] = offset_counts.sum(axis=1)
n_offset = offset_counts["Offset bifurcation types"].value_counts().sort_index()

onset_counts = metadata.groupby("patient_id")["bifurcation_start"].apply(
    lambda x: (
        x.value_counts()
        if len(x.value_counts()) > 0
        else pd.Series(name="count", index=["Unclear"], data=[None])
    )
)
onset_counts = onset_counts.unstack()
onset_counts["Onset bifurcation types"] = (onset_counts > 0).sum(axis=1)
onset_counts["Clear onset bifurcations"] = onset_counts.sum(axis=1)
n_onset = onset_counts["Onset bifurcation types"].value_counts().sort_index()

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(15, 6), sharey="row")

sns.barplot(
    onset_counts,
    x="Onset bifurcation types",
    y="Clear onset bifurcations",
    ax=axes[0],
    alpha=0.5,
)

for p, samples in zip(axes[0].patches, n_onset):
    axes[0].annotate(
        f"n={samples}",
        xy=(p.get_x() + p.get_width() / 2.0, 0),
        xytext=(0, 5),  # Offset text position above the bar
        textcoords="offset points",
        ha="center",
        va="bottom",
        fontsize=22,
    )

sns.barplot(
    offset_counts,
    x="Offset bifurcation types",
    y="Clear offset bifurcations",
    ax=axes[1],
    alpha=0.5,
)

for p, samples in zip(axes[1].patches, n_offset):
    axes[1].annotate(
        f"n={samples}",
        xy=(p.get_x() + p.get_width() / 2.0, 0),
        xytext=(0, 5),  # Offset text position above the bar
        textcoords="offset points",
        ha="center",
        va="bottom",
        fontsize=22,
    )

fig.savefig(FIG_PATH / "dynamical_types_vs_seizure_number.svg", dpi=300)
