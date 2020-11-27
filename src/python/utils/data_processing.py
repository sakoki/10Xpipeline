import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

def log_normalize(data, scaling_factor):
	"""Normalize gene expression data

	1. Dividing individual gene reads by total read count
	2. Multiply by scaling factor
	3. Take the natural log
	"""
	return np.log1p(((data.T) / data.T.sum()) * scaling_factor).T

def z_score(data, **kwargs):
	"""Z-score normalized data using sklearn StandardScaler"""
	scaler = StandardScaler(**kwargs)
	scaler.fit(data)
	scaled = scaler.transform(data)
	scaled = pd.DataFrame(scaled)
	scaled.columns = data.columns
	scaled.index = data.index
	return scaled