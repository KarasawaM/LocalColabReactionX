import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def make_2Dplot_from_df(df, df_xlabel, df_ylabel, plot_xlabel, plot_ylabel, outpath):
    """
    Example:
    make_plot(df=df,
            df_xlabel="index",
            df_ylabel="Delta E [kcal/mol]",
            plot_xlabel="Index",
            plot_ylabel=r"$\Delta E$ [kcal/mol]",
            outpath="plot.pdf")    
    """
    # Plot energies [kcal/mol]. read from df
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(df[df_xlabel], df[df_ylabel], marker='o', linestyle='-', linewidth=1.1, markersize=3)
    ax.set_xlabel(plot_xlabel, fontsize=13)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylabel(plot_ylabel, fontsize=13)
    ax.tick_params(axis='both', which='major', labelsize=12)
    plt.tight_layout()
    plt.savefig(outpath, bbox_inches='tight')

    return None


