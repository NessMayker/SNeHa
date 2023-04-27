"""
Normalizes an array.
"""
import numpy as np

def norm(probs):
    prob_factor = 1 / sum(probs)
    return [prob_factor * p for p in probs]

def normMap(data, percentile_lowest = 1.5, percentile_highest = 98.5):
    p_low  = np.percentile(data, percentile_lowest)
    p_high = np.percentile(data, percentile_highest)
    normalized = np.clip(data, p_low, p_high)
    normalized = normalized - np.min(normalized)
    normalized /= np.max(normalized)
    return(normalized)
