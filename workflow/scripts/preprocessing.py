# -*- coding: utf-8 -*-
"""Retrive data from the UK Gov Website."""
import pandas as pd
from helpers.get_logger import get_logger

LOG = snakemake.log[0]  # noqa: F821
PARAMS = snakemake.params  # noqa: F821
OUTPUT = snakemake.output  # noqa: F821

logger = get_logger(__name__, LOG, True)

# Read Data
tests = pd.read_csv(
    PARAMS["tests"],
    usecols=["date", "newVirusTestsByPublishDate"],
    parse_dates=["date"],
    index_col="date",
)
logger.info(f"Data read from {PARAMS['tests']}")
logger.info(f"Tests: \n{tests.info()}")

cases = pd.read_csv(
    PARAMS["cases"],
    usecols=["date", "newCasesByPublishDate"],
    parse_dates=["date"],
    index_col="date",
)
logger.info(f"Data read from {PARAMS['cases']}")
logger.info(f"Cases: \n{cases.info()}")

deaths = pd.read_csv(
    PARAMS["deaths"],
    usecols=["date", "newDeaths28DaysByDeathDate"],
    parse_dates=["date"],
    index_col="date",
)
logger.info(f"Data read from {PARAMS['deaths']}")
logger.info(f"Deaths: \n{deaths.info()}")

# process data
data = (
    tests.join([cases, deaths], how="inner").assign(
        rollingAvgNewTests=lambda df: df.loc[:, ["newVirusTestsByPublishDate"]]
        .rolling(7)
        .mean(),
        rollingAvgNewCases=lambda df: df.loc[:, ["newCasesByPublishDate"]]
        .rolling(7)
        .mean(),
        newCasesPerTest=lambda df: df.newCasesByPublishDate
        / df.newVirusTestsByPublishDate,
        rollingAvgCasesPerTest=lambda df: df.loc[:, ["newCasesPerTest"]]
        .rolling(7)
        .mean(),
        rollingAvgNewDeaths=lambda df: df.loc[:, ["newDeaths28DaysByDeathDate"]]
        .rolling(7)
        .mean(),
        newDeathsPerCase=lambda df: df.newDeaths28DaysByDeathDate
        / df.newCasesByPublishDate,
        rollingAvgDeathsPerCase=lambda df: df.loc[:, ["newDeathsPerCase"]]
        .rolling(7)
        .mean(),
        newDeathsPerTest=lambda df: df.newDeaths28DaysByDeathDate
        / df.newVirusTestsByPublishDate,
        rollingAvgDeathsPerTest=lambda df: df.loc[:, ["newDeathsPerTest"]]
        .rolling(7)
        .mean(),
    )
    # Move date from index to column to ease plotting
    .reset_index()
)
logger.info(f"Processed data:\n{data.info()}")

data.to_csv(OUTPUT["data"])
