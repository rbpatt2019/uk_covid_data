# -*- coding: utf-8 -*-
"""A short script for visualising UK Covid data.

The URLs to download the data are provide,
and retrieve the most up-to-date information from
coronavirus.data.gov.uk.
"""
import matplotlib.pyplot as plt
import pandas as pd
from helpers.get_logger import get_logger

LOG = snakemake.log[0]  # noqa: F821
INPUT = snakemake.input  # noqa: F821
OUTPUT = snakemake.output  # noqa: F821

logger = get_logger(__name__, LOG, True)
plt.rcParams["figure.figsize"] = (18, 10)

# Read data
data = pd.read_csv(INPUT["data"], index_col=0, parse_dates=["date"])

# Plot case data
cases_mosaic = [["tests", "cases"], ["testPerCase", "testPerCase"]]
(_, cases_axes) = plt.subplot_mosaic(cases_mosaic, sharex=False, sharey=False)

_ = data.plot(
    x="date",
    y="newVirusTestsByPublishDate",
    kind="scatter",
    color="gray",
    s=8,
    ax=cases_axes["tests"],
    title="Virus tests conducted",
    xlabel="",
)
_ = data.plot(
    x="date",
    y="rollingAvgNewTests",
    kind="line",
    color="blue",
    legend=None,
    ax=cases_axes["tests"],
    xlabel="",
)

_ = data.plot(
    x="date",
    y="newCasesByPublishDate",
    kind="scatter",
    color="gray",
    s=8,
    ax=cases_axes["cases"],
    title="Cases by date reported",
    xlabel="",
)
_ = data.plot(
    x="date",
    y="rollingAvgNewCases",
    kind="line",
    color="blue",
    legend=None,
    ax=cases_axes["cases"],
    xlabel="",
)

_ = data.plot(
    x="date",
    y="newCasesPerTest",
    kind="scatter",
    color="gray",
    s=8,
    ax=cases_axes["testPerCase"],
    title="Test-rate adjusted new cases",
)
_ = data.plot(
    x="date",
    y="rollingAvgCasesPerTest",
    kind="line",
    color="blue",
    legend=None,
    ax=cases_axes["testPerCase"],
)

for ax in cases_axes.values():
    # First Omicron case detected in UK
    ax.axvline(pd.to_datetime("2021-11-27"), color="r", linestyle="--")
    # Lateral flow tests approved for home use (approximate)
    ax.axvline(pd.to_datetime("2020-12-20"), color="g", linestyle="--")

plt.savefig(OUTPUT["cases"], dpi=300, bbox_inches="tight")
plt.close()

# Plot Deaths
# Use gridspec because layout is non-columnar
fig = plt.figure(constrained_layout=True)
gs = fig.add_gridspec(2, 6, wspace=0.5)

# Create axes
ax1 = fig.add_subplot(gs[0, :2])
ax2 = fig.add_subplot(gs[0, 2:4])
ax3 = fig.add_subplot(gs[0, 4:])
ax4 = fig.add_subplot(gs[1, :3])
ax5 = fig.add_subplot(gs[1, 3:])

_ = data.plot(
    x="date",
    y="newVirusTestsByPublishDate",
    kind="scatter",
    color="gray",
    s=8,
    title="Virus tests conducted",
    xlabel="",
    ax=ax1,
)
_ = data.plot(
    x="date",
    y="rollingAvgNewTests",
    kind="line",
    color="blue",
    legend=None,
    xlabel="",
    ax=ax1,
)

_ = data.plot(
    x="date",
    y="newDeaths28DaysByDeathDate",
    kind="scatter",
    color="gray",
    s=8,
    title="Deaths within 28 days of positive test by date of death",
    xlabel="",
    ax=ax2,
)
_ = data.plot(
    x="date",
    y="rollingAvgNewDeaths",
    kind="line",
    color="blue",
    legend=None,
    xlabel="",
    ax=ax2,
)

_ = data.plot(
    x="date",
    y="newCasesByPublishDate",
    kind="scatter",
    color="gray",
    s=8,
    title="Cases by date reported",
    xlabel="",
    ax=ax3,
)
_ = data.plot(
    x="date",
    y="rollingAvgNewCases",
    kind="line",
    color="blue",
    legend=None,
    xlabel="",
    ax=ax3,
)

_ = data.plot(
    x="date",
    y="newDeathsPerTest",
    kind="scatter",
    color="gray",
    s=8,
    title="Test-rate adjusted new deaths",
    ax=ax4,
)
_ = data.plot(
    x="date",
    y="rollingAvgDeathsPerTest",
    kind="line",
    color="blue",
    legend=None,
    ax=ax4,
)

_ = data.plot(
    x="date",
    y="newDeathsPerCase",
    kind="scatter",
    color="gray",
    s=8,
    title="Case-rate adjusted new deaths",
    ax=ax5,
)
_ = data.plot(
    x="date",
    y="rollingAvgDeathsPerCase",
    kind="line",
    color="blue",
    legend=None,
    ax=ax5,
)


for ax in [ax1, ax2, ax3, ax4, ax5]:
    # First Omicron case detected in UK
    ax.axvline(pd.to_datetime("2021-11-27"), color="r", linestyle="--")
    # Lateral flow tests approved for home use (approximate)
    ax.axvline(pd.to_datetime("2020-12-20"), color="g", linestyle="--")

plt.savefig(OUTPUT["deaths"], dpi=300, bbox_inches="tight")
plt.close()
