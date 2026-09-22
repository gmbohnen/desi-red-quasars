import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import gaussian_kde


def contour_scatter_plot(df,x,y="z-w3",y_threshold=3.9,y_threshold_label="ERQ threshold",subsample_col="rew_civ",subsample_cutoff=100,contour_levels=[25,55,85],figsize=(8,6),x_label=None,y_label=None,lower_x_lim=None,upper_x_lim=None,save_path=None,x_log=False,y_log=False,ax=None,subsample_color="red",threshold_color="DodgerBlue",show=False):
    density_df = df[[x,y]].copy()

    coords = density_df.values.T
    kde = gaussian_kde(coords)
    density_df["density"] = kde(coords)

    scatter_cutoff = np.percentile(density_df["density"],contour_levels[0])

    scatter_part = density_df[density_df["density"] < scatter_cutoff]
    # contour_part = density_df[density_df["density"] >= scatter_cutoff]

    x_vals = density_df[x].values
    y_vals = density_df[y].values

    xx, yy = np.mgrid[x_vals.min():x_vals.max():200j, y_vals.min():y_vals.max():200j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    grid_density = kde(positions).reshape(xx.shape)

    ## PLOTTING
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    ax.contour(xx, yy, grid_density, levels=[np.percentile(density_df["density"], val) for val in contour_levels], colors="black", linewidths=0.8)
    sns.scatterplot(data=scatter_part,x=x,y=y,s=1.5,alpha=0.9,ax=ax,color="black")

    if subsample_col is not None:
        temp = df[df[subsample_col] > subsample_cutoff]
        sns.scatterplot(data=temp,x=x,y=y,s=4,alpha=0.9,ax=ax,color=subsample_color)

    if y_threshold is not None:
        ax.axhline(y_threshold,label=y_threshold_label,linestyle="--",color=threshold_color)

    ax.set_xlabel(x_label if x_label is not None else x)
    ax.set_ylabel(y_label if y_label is not None else y)

    if lower_x_lim is None:
        lower_x_lim, _ = ax.get_xlim()
    if upper_x_lim is None:
        _, upper_x_lim = ax.get_xlim()
    
    ax.set_xlim(lower_x_lim,upper_x_lim)

    if x_log:
        ax.set_xscale("log")
    if y_log:
        ax.set_yscale("log")

    if ax is None:
        plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path)
    
    if show:
        plt.show()
