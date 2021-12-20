rule plotting:
    input:
        data=rules.preprocessing.output.data,
    output:
        cases=report(
            "results/plotting/adjusted_cases.png",
            caption="../reports/adjusted_cases.rst",
            category="Adjusted Cases",
        ),
        deaths=report(
            "results/plotting/adjusted_deaths.png",
            caption="../reports/adjusted_deaths.rst",
            category="Adjusted Deaths",
        ),
    log:
        "results/logs/plotting/plotting.log",
    benchmark:
        "results/benchmarks/plotting/plotting.txt"
    conda:
        "../envs/plotting.yaml"
    script:
        "../scripts/plotting.py"
