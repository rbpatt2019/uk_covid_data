from snakemake.utils import validate
from pathlib import Path
from yaml import safe_load

# Validate structures
validate(config, schema="../schema/config.schema.yaml")
