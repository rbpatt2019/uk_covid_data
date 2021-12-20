# Technical Documentation

What follows is a technical documentation of each step in the pipeline.

## Rule: preprocessing

Relatively straightforward.
The data are read in from the websites specified in the
[config file](config/config.yaml).
It is worth noting that certain columns are expected,
so any API change to the website may break the pipeline.
Once read in,
the desired columns are kept,
and the data joined together.
Test-rate and case-rate adjusted statistics are calculated
by dividing the apropriate column by the number of new tests.
To match the plot format from the
[government website][uk_covid],
7 day rolling averages are calculated for all statistics.

## Rule: plotting

The releveant plots are generated.
No data manipulations occurs.

[uk_covid]: https://coronavirus.data.gov.uk "UK Coronavirus Data"
