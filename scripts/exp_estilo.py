"""exp_estilo.py - EXPERIMENT, not part of the article: one style for every figure of the experiment.
Colours carry one meaning each, and the meaning is the one the reader already has:
  temperature   warm vermillion       precipitation  blue
  concentrates  dark red              diverts        dark blue      (normalized current, as in Fig. 2)
  both          amber: a bridge for one guild and a barrier for the other
  barrier       dark end of a sequential scale       open           light end
The pass ranks keep the colours of Fig. 1, and no other figure uses them. No grey bands mark zones: the figures show
where the barriers are, instead of saying it beforehand.
"""
import matplotlib
import matplotlib.patheffects
import matplotlib.pyplot as plt

GROUP = {"temperature": "#D55E00", "precipitation": "#0072B2"}
RANK = {"continental": "#7B3294", "terminal": "#008837", "tributary": "#8C8C8C"}   # Fig. 1 ranks, no colour reused
CONC, DIV = "#B2182B", "#2166AC"                                                  # normalized current, as in Fig. 2
BOTH = "#E69F00"                               # a pass that is a bridge for one guild and a barrier for the other
BARRIER = "magma"                                                                 # zone conductance: dark is a barrier
HALO = [matplotlib.patheffects.withStroke(linewidth=1.6, foreground="w")]
LAT_TICKS = range(-20, 11, 5)
VLAB = {"bio01": "BIO1 mean temperature", "bio02": "BIO2 diurnal range", "bio03": "BIO3 isothermality",
        "bio04": "BIO4 temperature seasonality", "bio05": "BIO5 max. temp. of warmest month",
        "bio06": "BIO6 min. temp. of coldest month", "bio07": "BIO7 temperature annual range",
        "bio08": "BIO8 temp. of wettest quarter", "bio09": "BIO9 temp. of driest quarter",
        "bio10": "BIO10 temp. of warmest quarter", "bio11": "BIO11 temp. of coldest quarter",
        "bio12": "BIO12 annual precipitation", "bio13": "BIO13 precip. of wettest month",
        "bio14": "BIO14 precip. of driest month", "bio15": "BIO15 precipitation seasonality",
        "bio16": "BIO16 precip. of wettest quarter", "bio17": "BIO17 precip. of driest quarter",
        "bio18": "BIO18 precip. of warmest quarter", "bio19": "BIO19 precip. of coldest quarter"}


VSHORT = {"bio01": "BIO1 mean temp.", "bio02": "BIO2 diurnal range", "bio03": "BIO3 isothermality",
          "bio04": "BIO4 temp. seasonality", "bio05": "BIO5 warmest month", "bio06": "BIO6 coldest month",
          "bio07": "BIO7 annual range", "bio08": "BIO8 wettest quarter", "bio09": "BIO9 driest quarter",
          "bio10": "BIO10 warmest quarter", "bio11": "BIO11 coldest quarter", "bio12": "BIO12 annual precip.",
          "bio13": "BIO13 wettest month", "bio14": "BIO14 driest month", "bio15": "BIO15 precip. seasonality",
          "bio16": "BIO16 wettest quarter", "bio17": "BIO17 driest quarter", "bio18": "BIO18 warmest quarter",
          "bio19": "BIO19 coldest quarter"}


def setup():
    plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.5, "xtick.major.width": 0.5,
                         "ytick.major.width": 0.5, "pdf.fonttype": 42, "savefig.dpi": 600})


def lat_label(v):
    return f"{abs(v)}°{'S' if v < 0 else 'N' if v else ''}"


def axis_lat(ax, lo, hi, labels=True):
    ax.set_ylim(lo, hi); ax.set_yticks(LAT_TICKS)
    ax.set_yticklabels([lat_label(v) for v in LAT_TICKS] if labels else [])
    ax.spines[["top", "right"]].set_visible(False)


def panel(ax, letter, text=""):
    ax.annotate(letter, (0, 1), xycoords="axes fraction", xytext=(0, 3), textcoords="offset points", fontsize=8,
                fontweight="bold", va="bottom", ha="left")
    if text:
        ax.annotate(text, (0, 1), xycoords="axes fraction", xytext=(9, 3), textcoords="offset points", fontsize=6.5,
                    va="bottom", ha="left")


def save(fig, path):
    fig.savefig(path + ".png", bbox_inches="tight", pad_inches=0.03)
    fig.savefig(path + ".pdf", bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("->", path)
