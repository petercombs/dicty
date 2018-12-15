import matplotlib.pyplot as mpl
import pandas as pd
from scipy import stats
from numpy import log10

if __name__ == "__main__":
    stalk = pd.read_table(
        "analysis/results/combined.Stalk.closest_mutant.bed",
        header=None,
        na_values=[".", "-1"],
    )
    spore = pd.read_table(
        "analysis/results/combined.Spore.closest_mutant.bed",
        header=None,
        na_values=[".", "-1"],
    )
    random = pd.read_table(
        "analysis/results/combined.Random.closest_mutant.bed",
        header=None,
        na_values=[".", "-1"],
    )
    data = [
        list(stalk.iloc[:, 9].dropna().apply(log10)),
        list(spore.iloc[:, 9].dropna().apply(log10)),
        list(random.iloc[:, 9].dropna().apply(log10)),
    ]
    mpl.violinplot(data, showextrema=False, showmedians=True)
    stalk_test = stats.mannwhitneyu(data[0], data[-1], alternative="less")
    spore_test = stats.mannwhitneyu(data[1], data[-1], alternative="less")
    stalk_sig = "*" if stalk_test.pvalue < .05 else ""
    spore_sig = "*" if spore_test.pvalue < .05 else ""
    print("Spore: ", spore_test)
    print("Stalk: ", stalk_test)
    mpl.xticks([1, 2, 3], ["Stalk" + stalk_sig, "Spore" + spore_sig, "Random"])
    mpl.ylabel("Log10 Distance\nto Nearest Mutant")
    mpl.tight_layout()
    mpl.savefig("analysis/results/mutant_distance.png")
