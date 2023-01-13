import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt



def calculate_dn_ds(n, s, N, S):
    pn = float(n) / float(N)
    ps = float(s) / float(S)
    #dn = -3.0 / 4.0 * (np.log(1 - ((4.0 / 3.0) * pn)))
    #ds = -3.0 / 4.0 * (np.log(1 - ((4.0 / 3.0) * ps)))
    dn = np.log(1 + pn)
    ds = np.log(1 + ps)
    return dn / ds


def plot_dn_ds(data_to_plot, factor, filename, title, legend=True, yrange=None):
    fig = plt.figure(figsize=(20, 8))
    sns.lineplot(
        data=data_to_plot, 
        x="month", 
        y="dn_ds", 
        hue=factor, 
        legend=False, 
        linestyle="--", 
        size=1,
        palette='colorblind')
    sns.scatterplot(
        data=data_to_plot, 
        x="month", 
        y="dn_ds", 
        hue=factor,
        s=100,
        legend=legend,
        style=factor, 
        palette='colorblind')
    plt.xticks(rotation=30)
    plt.title(title)
    plt.ylabel("dN/dS")
    plt.xlabel(None)
    plt.grid(axis="y")
    if yrange:
        plt.ylim(yrange)
    if legend:
        ax = fig.axes[0]
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(
            #handles=handles[1:], 
            #labels=labels[1:],
            loc='upper left'
        )
    sns.despine(left=True)
    plt.savefig(filename)
    
    
def plot_dn_ds_overall(data, N, S, title, filename):
    
    data_to_plot = data.groupby(["month", "source"]).sum().reset_index().sort_values("month")
    data_to_plot["dn_ds"] = data_to_plot[["s", "ns"]].apply(
        lambda x: calculate_dn_ds(s=x[0], n=x[1], S=S, N=N), axis=1)
    
    plot_dn_ds(
        data_to_plot=data_to_plot, 
        factor="source", 
        filename=filename, 
        legend=True,
        title=title
    )
    

def plot_dn_ds_by_country(data, title, filename, N, S, legend=True):
    
    data_to_plot = data.groupby(["month", "country"]).sum().reset_index().sort_values("month")
    data_to_plot["dn_ds"] = data_to_plot[["s", "ns"]].apply(
        lambda x: calculate_dn_ds(s=x[0], n=x[1], S=S, N=N), axis=1)
    
    plot_dn_ds(
        data_to_plot=data_to_plot, 
        factor="country", 
        filename=filename, 
        legend=legend,
        title=title
    )
    
    
def plot_dn_ds_by_gene(data, data_fractions, title, filename, legend=True, yrange=None):
    
    data_to_plot = data.groupby(["month", "region_name"]).sum().reset_index().sort_values("month")
    
    data_to_plot["dn_ds"] = data_to_plot[["s", "ns", "region_name"]].apply(
        lambda x: calculate_dn_ds(
            s=x[0], 
            n=x[1], 
            S=data_fractions[data_fractions.gene == x[2]].S.iloc[0], 
            N=data_fractions[data_fractions.gene == x[2]].NS.iloc[0]), axis=1)
    
    plot_dn_ds(
        data_to_plot=data_to_plot, 
        factor="region_name", 
        filename=filename, 
        legend=legend, 
        yrange=yrange,
        title=title
    )
    
    
def plot_dn_ds_by_domain(data, data_fractions, title, filename, legend=True, yrange=None):
    
    data_to_plot = data.groupby(["month", "region_name"]).sum().reset_index().sort_values("month")
    
    data_to_plot["dn_ds"] = data_to_plot[["s", "ns", "region_name"]].apply(
        lambda x: calculate_dn_ds(
            s=x[0], 
            n=x[1], 
            S=data_fractions[data_fractions.domain == x[2]].S.iloc[0], 
            N=data_fractions[data_fractions.domain == x[2]].NS.iloc[0]), axis=1)
    
    plot_dn_ds(
        data_to_plot=data_to_plot, 
        factor="region_name", 
        filename=filename, 
        legend=legend, 
        yrange=yrange,
        title=title
    )