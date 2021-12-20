# uk_covid_data - A Data Analysis Snakemake Workflow

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CI/CD](https://github.com/rbpatt2019/uk_covid_data/actions/workflows/cicd.yaml/badge.svg)](https://github.com/rbpatt2019/uk_covid_data/actions/workflows/cicd.yaml)
[![Codestyle: Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Codestyle: snakefmt](https://img.shields.io/badge/code%20style-snakefmt-000000.svg)](https://github.com/snakemake/snakefmt)

A [Snakemake][sm] workflow for visualising UK Covid-19 data.

> This work was done during my spare time.
> Any opinions expressed herein are exclusively **_mine_**.
> The aforementioned opinions do **_not_** represent the opinions
> of any employer or organisation with which I may be associated.

## Motivation

While reading about the current Covid-19/Omicron situation in the UK,
I realised that nearly all,
if not all,
of the data reported was "raw counts".
That is,
it did not consider the rate of testing.
Given how strongly the rate of testing will impact the number of cases reported -
e.g. a 10x increase in testing would (likely) result in a 10x increases in cases,
all else being equal -
this surprised me!
So I decided to make my own plots!
This pipeline represents as unbiased an attempt as possible
to investigate the data reported on the UK's
[Coronavirus website][uk_covid].

To be consistent with the Government's dashboard,
each plot shows the counts for the day,
as well as a 7-day rolling average.
Additionally,
I've added vertical lines for dates of interest -
namely,
the approval of the use of lateral flow tests
and the first detection of the omicron variant.

All the usual data science best practices still apply!
I've pipelined it for reproduciblility
and setup CI/CD to maintain standards.
Finally,
by taking advantage of the integrated [Conda][conda] and [Singularity][sing] support,
we can run the whole thing in an isolated environment.

## Notes on Installation

### Necessary Software

This pipeline needs [conda][conda]
and [snakemake][sm]
installed,
and runs best if you also have [singularity][sing]
installed,
though it's not required.

Snakemake recommends using [mambaforge][mambaforge]
as your base conda,
which I would also recommend.
Installation instructions are at the above link.
If you prefer a vanilla conda installation,
you can always try `mamba` following the instructions at the above snakemake link.
Once you have conda installed,
install snakemake as outlined on their page
(again, see the above link)
and activate your snakemake environment.

If you are running on a Linux system,
then singularity can be installed from conda like so:

```shell
conda install -n snakemake -c conda-forge singularity
```

It's a bit more challenging for other operating systems.
Your best bet is to follow their instructions
[here][sing_install].
But don't worry!
**Singularity is _not_ regquired!**
Snakemake will still run each step in its own Conda environment,
it just won't put each Conda environment in a container.

### Get the Source Code

Navigate to the [release][releases]
page on github and download the most recent version.
The following will do the trick:

```shell
curl -s https://api.github.com/repos/rbpatt2019/uk_covid_data/releases/latest |
grep tarball_url |
cut -d " " -f 4 |
tr -d '",' |
xargs -n1 curl -sL |
tar xzf -
```

After querying the github api to get the most recent release information,
we grep for the desired URL,
split the line and extract the field,
trim superfluous characters,
use `xargs` to pipe this to `curl` while allowing for re-directs,
and un-tar the files.
Easy!

Alternatively,
for the bleeding edge,
please clone the repo like so:

```shell
git clone https://github.com/rbpatt2019/uk_covid_data
```

> :warning: **Heads Up!**
> The bleeding edge may not be stable,
> as it contains all active development.

However you choose to install it,
`cd` into the directory.

### Running

Once you've installed the above software,
and fetched the code,
running the pipeline is as simple as:

```shell
snakemake --use-conda --use-singularity --cores 6
```

If you aren't using `singularity`,
then leave off the apropriate flag, as so:

```shell
snakemake --use-conda --cores 6
```

And `snakemake` will automatically leave it off.

Once the pipeline is completed,
run:

```shell
snakemake --report results/report.html
```

to generate a documented html file displaying all the results
and relevant runtime information.

## Notes on Configuration

To change the number of cores Snakemake uses,
change the number passed to the `--cores` flag,
like so:

```shell
snakemake --use-conda --use-singularity --cores 2
```

Snakemake gets cranky if that flag is omitted,
however,
so don't forget it!

> :warning:  **Be sure to change the configuration to suit your project!**

For a full discussion of configuration,
please see the [configuration README](config/README.md).

Briefly,
the general configuration file must be located at `config/config.yaml`.
This is schema validated,
so invalid entries will block runs.

## Notes on Data

The data used here is pulled from the UK's
[Coronavirus website][uk_covid].
The data will be pulled fresh when the pipeline is run,
guaranteeing up-to-date results.

For testing data,
total have been chosen over lab tests,
as a significant fraction -
somewhere around 50%, based on a quick glance at the website -
are home-performed tests.

Deaths are counted as
"deaths within 28 days of positive test by date of death"
while cases are counted as
"cases by date reported".
Though there are other options,
I feel these two are most likely to be truly representative.

## Notes on the tools

The analysis pipeline was run using Snakemake v6.12.3.
The full version and software lists can be found under the relevant yaml files in `workflow/envs`.
The all reasonable efforts have been made to ensure that the repository adheres to the best practices
outlined [here](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html).

## Notes on the analysis

For a full discussion on the analysis methods,
please see the [technical documentation](workflow/documentation.md).

Briefly,
data were pulled from the [government website][uk_covid],
followed by brief processing to join all data.
Then,
the necessary test-corrected and rolling metrics were calculated.
Finally, the data were plotted.

[sm]: https://snakemake.readthedocs.io/en/stable/index.html "Snakemake"
[uk_covid]: https://coronavirus.data.gov.uk "UK Coronavirus Data"
[conda]: https://docs.conda.io/en/latest/ "Conda"
[sing]: https://sylabs.io/singularity/ "Singularity"
[mambaforge]: https://github.com/conda-forge/miniforge#mambaforge "Mambaforge"
[sing_install]: https://sylabs.io/guides/3.8/admin-guide/installation.html#installation-on-windows-or-mac "Singularity Install"
[releases]: https://github.com/IMS-Bio2Core-Facility/single_snake_sequencing/releases "Releases"
