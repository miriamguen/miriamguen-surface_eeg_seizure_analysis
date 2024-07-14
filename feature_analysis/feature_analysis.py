"This code will be used for the analysis of the features extracted from the components"
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

#%% define data version and paths
DATA_ALIAS = 'example_data'
# folder definition
BASE_PATH = Path(os.getcwd().split('feature_analysis')[0])
DATA_DIR = BASE_PATH / DATA_ALIAS /  'features'
FIG_DIR = DATA_DIR / 'figures'
os.makedirs(FIG_DIR, exist_ok=True)
#%% define plotting defaults
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
#%% Load the onset seizure features

file_list = list(filter(lambda x: x.endswith('.mat_start_bif_variance_both.csv') , os.listdir(DATA_DIR / 'seizure')))

seizure_features_list = []
for file in file_list:
    seizure_features_list.append(pd.read_csv(DATA_DIR / 'seizure'/ file))


seizure_features = pd.concat(seizure_features_list)
del seizure_features_list