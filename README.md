# Snakemake Pipeline for Impact Simulation

This Snakemake pipeline is designed to simulate and aggregate the impact of events using the Climada library. It takes input data related to direct impacts, supply chain information, and other configuration parameters to generate simulations and aggregate the results.

## Prerequisites

    Python >= 3.6
    Snakemake >= 6.0
    Climada library
    Conda environment specified in envs/climada_env.yml

## Getting Started

    Clone this repository to your local machine.
    Make sure you have the required dependencies installed.
    Modify the config/config.yaml file to configure the pipeline settings according to your needs.

## Pipeline Structure

The pipeline consists of several rules that perform specific tasks. Here's an overview of each rule:

### Rule: impact_csv_csr_to_pkl
This rule reads CSV and NPZ files produced with climada and containing direct impact data and creates an Impact object, which is then saved using pickle format.

### Rule: create_supchain
This rule creates a supply chain object, calculates direct impacts, and generates meta information about
events using the Climada library. The output includes a supply chain object, impact dataframes, and other relevant data.

### Rule: simulate
This rule runs simulations for individual events or sets of events based on the supply chain and impact data. The simulations are performed using the Climada library and the results are saved in separate Parquet files.

### Rule: reduce
This rule aggregates the results of isolated events simulations. It combines the results for each variable across all events and saves them as separated Parquet files.

### Rule: all_results
This rule collects all simulation results for further processing.

### Rule: meta_df
This rule aggregates and produces a meta dataframe that groups all simulation results. The meta dataframe provides a comprehensive overview of the impact simulation outcomes.

## Required Input Files

Before running the pipeline, ensure that you have the following input files prepared and placed in the appropriate directories:
1. Direct Impact Data

Direct impact data files consist of CSV and NPZ pairs representing the impacts of events. The CSV file contains event-specific information, and the NPZ file contains impact values.

CSV File Format:

    Path: inputs/direct_impacts_yearsets/sample_{n}_5y_impact_na_wp.csv
    Columns: EventID, EventDate, Region, Sector, ImpactValue, etc.

NPZ File Format:

    Path: inputs/direct_impacts_yearsets/sample_{n}_5y_impact_na_wp.npz
    Contents: Direct impact values as NumPy arrays

2. MRIOT Data

The MRIOT (Multi-Regional Input-Output Table) data is used to create a supply chain object and perform impact calculations. The MRIOT file should be in pickle format.

MRIOT File Format:

    Path: inputs/EXIOBASE3-2010/mriot_original.pkl

3. Sector Aggregation Data

The sector aggregation data is used to aggregate sectors within the supply chain object.

Sector Aggregation File Format:

    Path: inputs/sectors_common_aggreg.ods

4. Additional Configurations

The pipeline relies on various configuration parameters specified in the config/config.yaml file. Ensure that this file is correctly configured with the necessary values for parameters such as rebuild_tau, duration, order_type, alpha_tau, psi, monetary_factor, step_to_eq, and variables.
Folder Structure

To maintain a standardized folder structure, consider organizing your input files and outputs as follows:

    project_root/
    │
    ├─ config/
    │   └─ config.yaml
    │
    ├─ envs/
    │   └─ climada_env.yml
    │
    ├─ inputs/
    │   ├─ direct_impacts_yearsets/
    │   │   ├─ sample_{n}_5y_impact_na_wp.csv
    │   │   ├─ sample_{n}_5y_impact_na_wp.npz
    │   │   └─ ...
    │   │
    │   ├─ EXIOBASE3-2010/
    │   │   └─ mriot_original.pkl
    │   │
    │   └─ sectors_common_aggreg.ods
    │
    ├─ outputs/
    │   ├─ ...
    │   │
    │   └─ simulations/
    │       ├─ sample_{n}/
    │       │   ├─ EXIOBASE3-2010_supchain.pkl
    │       │   ├─ EXIOBASE3-2010_df_impact.parquet
    │       │   ├─ EXIOBASE3-2010_meta_df_impact.parquet
    │       │   ├─ ...
    │       │
    │       └─ ...
    │
    ├─ scripts/
    │   ├─ yearset_to_pkl.py
    │   ├─ supchain.py
    │   ├─ indirect_impact.py
    │   └─ meta_df_building.py
    │
    ├─ Snakefile
    └─ README.md

Make sure to replace {n} and {mriot_name} with actual values in your filenames and directory paths.

## Running the Pipeline

To run the pipeline, execute the following command in the terminal in a python environment with snakemake installed:

    snakemake meta_df --use-conda -p --rerun-incomplete

Replace <num_cores> with the number of cores you want to use for parallel processing.
Replace `meta_df` by any other rule (or output file) if you want to run just this part.
You can also add a `> "./snakelog-$(date +"%FT%H%M%z").log" 2>&1 &` part to create a 
timestamped log file of snakemake output and keep the terminal handle.
Note that you do not require the other dependencies to be in that environment, 
as snakemake will create the required environments by itself. Also note that this requires a bit of diskspace.

I strongly suggest you use a cluster/computation server. In this case, setup a snakemake profile in your
home folder `/home/{user}/.config/snakemake/{profile_name}/config.yaml`. 
I use the following setup (based from https://github.com/jdblischak/smk-simple-slurm) for a slurm based cluster, 
make sure you also get `status-sacct.sh` from this repo and add it in the same folder as `config.yaml`.
    
    cluster:
      mkdir -p /scratchu/sjuhel/logs/smk/{rule} &&
      sbatch
        --parsable
        --mem={resources.mem_mb}
        --job-name=smk-{rule}-{wildcards}
        --cpus-per-task={threads}
        --output=/scratchu/sjuhel/logs/smk/{rule}/{rule}-{wildcards}-%j.out
        --time={resources.time}
        --partition={resources.partition}
    default-resources:
      - mem_mb=2000
      - partition=zen16
      - time=60
      - threads=4
    restart-times: 0
    max-jobs-per-second: 10
    max-status-checks-per-second: 1
    local-cores: 2
    latency-wait: 60
    jobs: 16
    keep-going: False
    rerun-incomplete: True
    printshellcmds: True
    scheduler: greedy
    use-conda: True
    conda-prefix: {where you want to store the python envs created by snakemake}
    conda-base-path: {your python install path (/home/user/conda for example)}
    cluster-status: status-sacct.sh

You can then simply invoke snakemake with:

    snakemake meta_df --profile simple > "./snakelog-$(date +"%FT%H%M%z").log" 2>&1 &

## Output

The pipeline generates various output files in the outputs/ directory:

  1. Direct impact data
     - a climada Impact object in pkl format `outputs/sample_{n}/direct_impact.pkl`
     - a climada SupplyChain object in pkl format `outputs/sample_{n}/{mriot_name}_supchain.pkl`
     - the climada compatible mriot in pkl format (then used for simulation) `outputs/sample_{n}/{mriot_name}.pkl`
     - a dataframe (parquet) of the direct impacts with (region,sector) as column, events as row `outputs/sample_{n}/{mriot_name}_df_impact.parquet`.
     - a dataframe (parquet) of some aggregated infos and variables on each events (total direct damages, affected regions, recovery duration) `outputs/sample_{n}/{mriot_name}_meta_df_impact.parquet`.
  2. Supply chain objects
  3. Simulated results for each event and variable
     - These are stored in "outputs/simulations/sample_{n}/{mriot_name}_{inv_sce}_event_{event_id}/{variable}.parquet"
     - These parquets files are dataframes that contains the values of the variable at each step of the simulation.
       {event_id} is either the id of a specific event (then the file is the results for that isolated event) or
       "separated" (then the file is the aggregation of all isolated events), or "aggregated" (then the file is
       the results for all events considered together)
  4. Meta dataframe summarizing the simulation outcomes
     - Stored in "outputs/results_meta_df.parquet"
     - Shows the total direct and total indirect impact for each separated and aggregated simulation, for the different inventory parametrisation, the different variables, on a per region/sector basis

## Configuration

Adjust the pipeline configuration in the config/config.yaml file to specify input paths, parameter values, and other settings necessary for the simulation.
Acknowledgments

This pipeline was developed using the Snakemake workflow management system, the Climada platform and the BoARIO model.
References

    Snakemake Documentation: https://snakemake.readthedocs.io/
    Climada Documentation: https://climada-python.readthedocs.io/en/stable/
    BoARIO Documentation: https://spjuhel.github.io/BoARIO/
