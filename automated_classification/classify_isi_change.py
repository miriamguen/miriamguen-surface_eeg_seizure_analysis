import os
import pickle
from os.path import join
from pathlib import Path
import matplotlib as mpl
import pandas as pd
import shap
import matplotlib.pyplot as plt

from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

from classification_utils import *

mpl.use("svg")

EXPERIMENT_NAME = "detect_isi_change"

BASE_PATH = Path(os.getcwd())
FEATURE_ALIAS = "example_data"
DATA_PATH = BASE_PATH / FEATURE_ALIAS / "features"
RESULT_PATH = BASE_PATH / FEATURE_ALIAS / "classification_results"
FIG_PATH = BASE_PATH / "figures" / "classification"

RESULT_PATH.mkdir(parents=True, exist_ok=True)
FIG_PATH.mkdir(parents=True, exist_ok=True)


MODEL_PARAMS = {
    "max_depth": 3,
    "n_estimators": 500,
    "random_state": 12,
    "class_weights": {0: 0.1, 1: 0.8},  # use the class weights
    "n_jobs": -1,
}

# %% load the relevant data

onset_data = pd.read_csv(
    DATA_PATH / "all_onset_features_for_classification.csv", index_col="index"
).rename(columns={"bifurcation_start": "label"})
onset_data["onset_offset"] = 0

offset_data = pd.read_csv(
    DATA_PATH / "all_offset_features_for_classification.csv", index_col="index"
).rename(columns={"bifurcation_end": "label"})
offset_data["onset_offset"] = 1

raw_data = pd.concat([onset_data, offset_data], axis=0, join="inner").reset_index(
    drop=True
)
# drop all non clear
raw_data = raw_data[raw_data["label"] != "not_clear"]
print(raw_data["label"].value_counts())

raw_data = raw_data.assign(
    label=raw_data["label"].apply(lambda x: 1 if "SNIC" in x else 0)
)

print(
    f"out of {len(raw_data)} examples there are {raw_data['label'].sum()} isi slowing examples"
)


print("stop")
# pre_process data
data, metadata, processing_info = prepare_data(raw_data, None)


# %% run the leave one out cross validation
metadata, train_results_all, test_results_all, overall_test_performance = (
    loo_training_and_validation(data, metadata, MODEL_PARAMS)
)

# save artifacts

results = test_results_all.describe()
results.loc["overall", :] = overall_test_performance
os.makedirs(RESULT_PATH / EXPERIMENT_NAME, exist_ok=True)

test_results_all.to_csv(RESULT_PATH / EXPERIMENT_NAME / f"loo_results.csv")
results.to_csv(RESULT_PATH / EXPERIMENT_NAME / f"loo_results_summary.csv")
data.to_csv(RESULT_PATH / EXPERIMENT_NAME / "analysis_data.csv")
metadata.to_csv(RESULT_PATH / EXPERIMENT_NAME / "analysis_metadata_with_loo_pred.csv")


# %% Analyze model explainability
# train model on all data for explainability analysis, save the trained model
model = fit_model(data, metadata["label"], model_params=MODEL_PARAMS)

with open(
    RESULT_PATH / EXPERIMENT_NAME / f"identify_clear_seizure_transition_model.pkl", "wb"
) as file:
    pickle.dump(model, file)


clear_probability = lambda x: model.predict_proba(x)[:, 1]
med = pd.DataFrame(data.median(axis=0)).T

explainer = shap.Explainer(
    clear_probability,
    med,
    algorithm="permutation",
    seed=0,
    output_names=["Unclear", "Clear"],
)

shap_values = explainer(data)
shap_values_df = pd.DataFrame(
    data=shap_values.values, columns=data.columns, index=data.index
)
shap_values_df.to_csv(RESULT_PATH / EXPERIMENT_NAME / "shap_values_all.csv")

# %% plot results

fig1 = plt.figure()
shap.summary_plot(
    shap_values,
    data,
    feature_names=[x.replace("_", " ").capitalize() for x in data.columns],
    cmap="coolwarm",
    max_display=10,
)
fig1.suptitle("ISI slowing top features")
fig1.savefig(
    FIG_PATH / f"{EXPERIMENT_NAME}_summary_plot.svg",
    bbox_inches="tight",
    transparent=True,
)

fig2, ax = plt.subplots(1, 1)
preds = model.predict(data)

disp = ConfusionMatrixDisplay(
    confusion_matrix=confusion_matrix(
        metadata["label"], metadata["pred_loo"], normalize="true"
    ),
    display_labels=["No slowing", "Slowing"],
)

disp.plot(
    include_values=True,
    cmap="Blues",  # "coolwarm",
    ax=ax,
    xticks_rotation="horizontal",
    values_format=None,
    colorbar=False,
    im_kw={"alpha": 0.8},
)

# Define custom font properties for the matrix values
font_properties = {'fontsize': 16}

# Update the font properties of the text inside the confusion matrix
for text in disp.text_.ravel():
    text.set_fontsize(font_properties['fontsize'])



ax.set_yticklabels(["No slowing", "Slowing"], rotation=90)
ax.set_ylabel("Visual label", rotation=90, fontsize=20)
ax.set_xlabel("Classifier label", rotation=0, fontsize=20)

ax.set_xticklabels(["No slowing", "Slowing"], fontdict=font_properties)
ax.set_yticklabels(["No slowing", "Slowing"], rotation=90, fontdict=font_properties)

fig2.savefig(
    FIG_PATH / f"{EXPERIMENT_NAME}_confusion.svg",
    bbox_inches="tight",
    transparent=True,
)

print("all done")
