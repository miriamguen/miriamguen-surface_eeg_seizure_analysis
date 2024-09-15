import os
import pickle
from os.path import join
from pathlib import Path
import matplotlib as mpl
import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.preprocessing import RobustScaler
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_classif

from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from statsmodels.stats.inter_rater import cohens_kappa
from sklearn.metrics import confusion_matrix

mpl.use("svg")
# default parameters if not provided
MODEL_PARAMS = {
    "max_depth": 3,
    "n_estimators": 500,
    "random_state": 12,
    "class_weights": {0: 1, 1: 1},
    "n_jobs": -1,
}

# %% define help functions


def mean_bif_time(df):
    """replace unavailable bif times with the seizure detected mean

    Args:
        df (_type_): _description_

    Returns:
        _type_: _description_
    """
    mean_time = df.loc[df.bif_time != -1, "bif_time"].mean()
    df.loc[df.bif_time == -1, "bif_time"] = mean_time
    return df


def prepare_data(data, keep_n=None):
    """Prepare the feature data extracted for classification

    Args:
        data (_type_): the feature and labels
    """
    processing_info = {}
    data = data.assign(label=data.label.apply(int))
    data = data.assign(
        patient=data.file_name.apply(lambda x: "_".join(x.split("_")[0:3])).values
    )
    metadata = data[["patient", "label", "bif_time", "onset_offset"]]
    data["detacted_change_point"] = np.sign(data["bif_time"].values)
    data["bif_time"] = data.groupby("file_name", group_keys=False)[["bif_time"]].apply(
        mean_bif_time)

    data = data.drop(columns=["file_name", "comp_num", "patient", "label"])

    # remove features with over 5% nan values as they are unstable features for this case
    n = len(data) / 20
    keep_columns = data.isna().sum(axis=0).sort_values() < n

    data = data.loc[:, keep_columns]

    # # drop rows with over 5% nan values # The model shoud deal with nall values
    # drop_rows = data.isna().sum(axis=1) > (len(data.columns) / 20)

    # processing_info["drop_rows"] = drop_rows.sum()
    # print(f"Dropping {drop_rows.sum()} data rows due to noisy data")
    # data = data[~drop_rows].reset_index()
    # metadata = metadata[~drop_rows].reset_index()
    
    data = data.reset_index()
    metadata = metadata.reset_index()

    processing_info["feature_columns_names"] = data.columns
    print(
        f"feature matrix contains {len(data)} examples, out of them {metadata.label.sum()} are labeled as a clear transition"
    )

    processing_info["remaining_rows"] = len(metadata)
    processing_info["positive_labels"] = metadata.label.sum()
    processing_info["negative_labels"] = len(metadata) - metadata.label.sum()

    data = data.fillna(data.median())

    if "index" in data.columns:
        data = data.drop(columns=["index"])

    if keep_n is None:
        keep_columns = data.columns.values
    else:
        keep_n = len(data.columns)
        # select the top  keep_n features for a cleaner model
        selector = SelectKBest(f_classif, k=keep_n).fit(data, metadata.label)
        keep_columns = selector.get_feature_names_out()
        data = data[keep_columns]

    print(f"Keeping {len(keep_columns)} feature columns for task")
    processing_info["feature_columns"] = len(keep_columns)
    return (data, metadata, processing_info)


def fit_model(X_train, y_train, model_params):
    model = Pipeline(
        [
            (
                "scale",
                RobustScaler(),
            ),  # use robust scaler to for an outlier robust transformation
            (
                "rf",
                RandomForestClassifier(
                    max_depth=model_params["max_depth"],
                    n_estimators=model_params["n_estimators"],
                    random_state=model_params["random_state"],
                    class_weight=model_params["class_weights"],
                    n_jobs=model_params["n_jobs"],
                ),
            ),
        ]
    )

    return model.fit(X_train, y_train)


def produce_report(y_test, y_pred, y_probability):
    """_summary_

    Args:
        y_test (_type_): _description_
        y_pred (_type_): _description_
        y_probability (_type_): _description_

    Raises:
        ValueError: _description_

    Returns:
        _type_: _description_
    """
    report = {}

    cm = pd.DataFrame(confusion_matrix(y_test, y_pred, labels=[0, 1]))

    report["sample_number"] = len(y_pred)
    report["true_negative"] = cm.loc[0, 0]
    report["false_negative"] = cm.loc[1, 0]
    report["true_positive"] = cm.loc[1, 1]
    report["false_positive"] = cm.loc[0, 1]

    if (sum(y_test == 1) > 0) and (sum(y_test == 0) > 0):

        report["balanced_accuracy"] = balanced_accuracy_score(y_test, y_pred)
        report["roc_auc"] = roc_auc_score(y_test, y_probability)
        report["local_kappa"] = cohens_kappa(cm)["kappa"]
    elif (
        len(y_test == 1) > 0
    ):  # in cases with only clear seizure's no fp are possible, auc-roc not defined
        report["balanced_accuracy"] = np.mean(y_test == y_pred)
        report["roc_auc"] = np.nan
        report["local_kappa"] = cohens_kappa(cm)["kappa"]
    elif (
        len(y_test == 0) > 0
    ):  # in cases were no clear seizures exists roc auc is not defined
        report["balanced_accuracy"] = np.mean(y_test == y_pred)
        report["roc_auc"] = np.nan
        report["local_kappa"] = cohens_kappa(cm)["kappa"]
    else:
        raise ValueError("no data for this patient")
    return report


def eval_model(X_test, y_test, model):
    """_summary_

    Args:
        X_test (_type_): _description_
        y_test (_type_): _description_
        model (_type_): _description_

    Raises:
        ValueError: _description_

    Returns:
        _type_: _description_
    """
    y_pred = model.predict(X_test)
    y_probability = model.predict_proba(X_test)[:, 1]

    return produce_report(y_test, y_pred, y_probability)


def loo_training_and_validation(data, metadata, model_params=MODEL_PARAMS):
    """_summary_

    Args:
        data (_type_): _description_
        metadata (_type_): _description_

    Returns:
        _type_: _description_
    """
    logo = LeaveOneGroupOut()
    logo.get_n_splits(X=data, y=metadata["label"], groups=metadata["patient"])
    n_patients = len(set(metadata["patient"].values))
    train_results_all = {}
    test_results_all = {}

    for i, (train_index, test_index) in enumerate(
        logo.split(X=data, y=metadata["label"], groups=metadata["patient"])
    ):
        print(f"Round {i} - TRAIN: {len(train_index)} , TEST: {len(test_index)}")
        # divide the data into train and test for current split
        X_train, X_test = (
            data.loc[train_index, :],
            data.loc[test_index, :],
        )
        y_train, y_test = (
            metadata.loc[train_index, "label"],
            metadata.loc[test_index, "label"],
        )

        model = fit_model(X_train, y_train, model_params)
        train_results = eval_model(X_train, y_train, model)
        test_results = eval_model(X_test, y_test, model)
        # save run result
        train_results_all[i] = train_results
        test_results_all[i] = test_results
        # update test predictions in the metadata
        metadata.loc[test_index, "pred_prob_loo"] = model.predict_proba(X_test)[:, 1]
        metadata.loc[test_index, "pred_loo"] = model.predict(X_test)

    train_results_df = pd.DataFrame(train_results_all).T
    test_results_df = pd.DataFrame(test_results_all).T

    overall_test_performance = produce_report(
        metadata["label"], metadata["pred_loo"], metadata["pred_prob_loo"]
    )
    overall_test_performance["patient_nember"] = n_patients
    return (metadata, train_results_df, test_results_df, overall_test_performance)
