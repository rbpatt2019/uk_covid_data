# Configuration

The configuration keys that are expected are given below.
Don't worry about typos, etc.
These are all enforced with Snakemake's brilliant
[schema validation](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#validation).

## config.yaml

### preprocessing

- tests: `str, required`: URL from which to download the testing data.
- cases: `str, required`: URL from which to download the case data.
- deaths: `str, required`: URL from which to download the death data.
