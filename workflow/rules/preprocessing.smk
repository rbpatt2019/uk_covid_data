rule preprocessing:
    output:
        data="results/preprocessing/data.csv",
    params:
        tests=config["preprocessing"]["tests"],
        cases=config["preprocessing"]["cases"],
        deaths=config["preprocessing"]["deaths"],
    log:
        "results/logs/preprocessing/preprocessing.log",
    benchmark:
        "results/benchmarks/preprocessing/preprocessing.txt"
    conda:
        "../envs/preprocessing.yaml"
    script:
        "../scripts/preprocessing.py"
