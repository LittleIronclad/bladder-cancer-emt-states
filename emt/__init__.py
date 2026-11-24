"""EMT states in bladder cancer: NMF clustering, LASSO Cox survival modelling and validation."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

__version__ = "2.0.0"

plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 160, "savefig.bbox": "tight",
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
    "axes.spines.top": False, "axes.spines.right": False,
})


def savefig(fig, name):
    from .config import FIG_DIR
    fig.savefig(FIG_DIR / name)
    plt.close(fig)
    print(f"  figure → figures/{name}")


def banner(text):
    print("\n" + "=" * 70 + f"\n{text}\n" + "=" * 70)
