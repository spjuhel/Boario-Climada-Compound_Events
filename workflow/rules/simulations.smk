import pandas as pd


def get_all_parquets(wildcards):
    parquet = checkpoints.create_simulation_space.get(**wildcards).output
    simulations = pd.read_parquet(parquet)
    region_class = simulations[["mrio_region", "flood_class"]].drop_duplicates()
    tmp = "boario-global-floods/simulations/{parameters}/{mrio_region}~{flood_class}/parquets/{var}.parquet"
    tmp = expand(
        tmp, parameters=smk_parameters_space.instance_patterns, allow_missing=True
    )
    tmp = expand(tmp, var=config["variables_to_save"], allow_missing=True)
    tmp = expand(
        tmp,
        zip,
        mrio_region=region_class["mrio_region"],
        flood_class=region_class["flood_class"],
    )
    return tmp


rule all_parquets:
    input:
        get_all_parquets,


checkpoint create_simulation_space:
    input:
        rep_events=expand(
            "{flood_dir}/builded-data/{mrio_basename}/7_representative_events.parquet",
            flood_dir=config["flood_dir"],
            mrio_basename=smk_parameters_space.dataframe.mrio.str.split(
                "_", n=1, expand=True
            )[0].unique(),
        ),
        mrios=expand(
            "{mriot_dir}/pkls/{mrio_basename}/{mrio_basename}_{year}_full.pkl",
            mriot_dir=config["data-mriot"]["prefix"],
            mrio_basename=smk_parameters_space.dataframe.mrio.str.split(
                "_", n=1, expand=True
            )[0].unique(),
            year=smk_parameters_space.dataframe.mrio.str.split("_", expand=True)[
                1
            ].unique(),
        ),
    output:
        "boario-global-floods/simulation_space.parquet",
    benchmark:
        "benchmarks/simulation_space.csv"
    log:
        "logs/simulation_space.log",
    conda:
        "../envs/boario-global-floods.yml"
    resources:
        mem_mb=6000,
    script:
        "../scripts/simulation_space.py"


def get_simulation_inputs(wildcards):
    mrio_name = wildcards.mrio
    mrio_basename, year, agg = mrio_name.split("_")
    mrio_path = f"{config['data-mriot']['prefix']}/pkls/{mrio_basename}/{mrio_name}.pkl"
    sectors_scenario = f"{config['sectors_scenarios']['path']}/{mrio_basename}_{agg}_{wildcards.sectors_scenario}.csv"
    return {
        "mrio": mrio_path,
        "sectors_scenario": sectors_scenario,
    }


rule simulate:
    input:
        unpack(get_simulation_inputs),
        sim_df=rules.create_simulation_space.output,
    output:
        # format a wildcard pattern like "alpha~{alpha}/beta~{beta}/gamma~{gamma}"
        # into a file path, with alpha, beta, gamma being the columns of the data frame
        output_dir=directory(
            f"boario-global-floods/simulations/{smk_parameters_space.wildcard_pattern}/"
            + "{mrio_region}~{flood_class}"
        ),
        parquet_files=[
            f"boario-global-floods/simulations/{smk_parameters_space.wildcard_pattern}/"
            + "{mrio_region}~{flood_class}/parquets/"
            + var
            + ".parquet"
            for var in config["variables_to_save"]
        ],
        json_files=[
            f"boario-global-floods/simulations/{smk_parameters_space.wildcard_pattern}/"
            + "{mrio_region}~{flood_class}/jsons/"
            + json
            + ".json"
            for json in [
                "indexes",
                "equilibrium_checks",
                "simulated_events",
                "simulated_params",
            ]
        ],
        done=touch(
            f"boario-global-floods/simulations/{smk_parameters_space.wildcard_pattern}/"
            + "{mrio_region}~{flood_class}/sim.done"
        ),
    params:
        # automatically translate the wildcard values into an instance of the param space
        # in the form of a dict (here: {"alpha": ..., "beta": ..., "gamma": ...})
        sim_length=config["sim_length"],
    log:
        f"logs/simulations/{smk_parameters_space.wildcard_pattern}/"
        + "{mrio_region}~{flood_class}.log",
    threads: 4
    benchmark:
        (
            f"benchmarks/simulations/{smk_parameters_space.wildcard_pattern}/"
            + "{mrio_region}~{flood_class}.log"
        )
    resources:
        mem_mb=10000,
    conda:
        "../envs/boario-global-floods.yml"
    script:
        "../scripts/simulate.py"
