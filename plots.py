# plots.py
from __future__ import annotations
from typing import List, Dict
import pandas as pd
import plotly.express as px

def plot_gc_distribution(records: List[Dict]):
    df = pd.DataFrame(records)
    col = "GC_percent" if "GC_percent" in df.columns else ("GC_content" if "GC_content" in df.columns else None)
    if col is None:
        return None
    return px.histogram(df, x=col, nbins=30, title="GC% Distribution")

def plot_length_distribution(records: List[Dict]):
    df = pd.DataFrame(records)
    if "Length" not in df.columns:
        return None
    return px.histogram(df, x="Length", nbins=30, title="Sequence Length Distribution")

def plot_base_composition(records: List[Dict]):
    df = pd.DataFrame(records)
    needed = {"A","T","G","C","denom_used"}
    if not needed.issubset(df.columns):
        return None

    # mean fractions
    denom = df["denom_used"].replace({0: None})
    frac = pd.DataFrame({
        "A": df["A"] / denom,
        "T": df["T"] / denom,
        "G": df["G"] / denom,
        "C": df["C"] / denom,
    })
    means = frac.mean(skipna=True).reset_index()
    means.columns = ["Base", "MeanFraction"]
    return px.bar(means, x="Base", y="MeanFraction", title="Mean Base Fractions (A/T/G/C)", range_y=[0,1])

def plot_quality_distribution(records: List[Dict]):
    df = pd.DataFrame(records)
    if "Avg_quality" not in df.columns:
        return None
    return px.histogram(df, x="Avg_quality", nbins=30, title="Average Read Quality (Phred) Distribution")
