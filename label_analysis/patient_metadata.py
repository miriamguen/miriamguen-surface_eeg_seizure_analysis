"""Script used to summarize and print out the patient metadata"""

import pandas as pd
from pathlib import Path
import os


DATA_PATH = Path(os.getcwd()) / "example_data/source_data/metadata.csv"

data = pd.read_csv(DATA_PATH, index_col=0)

# Group by patient_id and aggregate
aggregated_data = (
    data.groupby("patient_id")
    .agg({"age": "first", "gender": "first", "lobe": "first", "etiology": "first"})
    .reset_index()
)

# Display the result
print(aggregated_data["age"].describe())
print(aggregated_data["gender"].value_counts())
print(aggregated_data["lobe"].value_counts())
print(aggregated_data["etiology"].value_counts())

print("done")
