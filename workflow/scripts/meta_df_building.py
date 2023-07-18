import sys
import logging
import pathlib
import pandas as pd

logger = logging.getLogger(__name__)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
fh = logging.FileHandler(snakemake.log[0])
fh.setLevel(logging.INFO)

def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return

    logger.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type, exc_value, exc_traceback),
            ]
        )
    )

sys.excepthook = handle_exception

logFormatter = logging.Formatter(
    "%(asctime)s [%(levelname)-5.5s] %(name)s %(message)s", datefmt="%H:%M:%S"
)
ch.setFormatter(logFormatter)
fh.setFormatter(logFormatter)
logger.addHandler(ch)
logger.addHandler(fh)

def process_impact_df(filepath, n, mriot):
    df = pd.read_parquet(filepath)
    df_sum = df.sum(axis=0).to_frame(name="total direct impact")
    df_sum['sample'] = n
    df_sum['MRIOT'] = mriot
    return df_sum

def process_res_df(filepath, n, mriot, simtype, var):
    df = pd.read_parquet(filepath)
    df_sum = (df-df.loc[0]).sum(axis=0).to_frame("total indirect impact")
    df_sum["sample"] = n
    df_sum["MRIOT"] = mriot
    df_sum["simtype"]= simtype
    df_sum["variable"] = var
    return df_sum

import re

mriots = []
ns = []

pattern_impact = r"sample_(\d+)/([^/]+)_df_impact.parquet"
pattern_results = r"sample_(\d+)/([^/]+)_event_(aggregated|separated)/([^/]+)\.parquet"
regex_impact = re.compile(pattern_impact)
regex_results = re.compile(pattern_results)

dfs_impact = []

for df_file in snakemake.input.impacts:
    match = regex_impact.search(df_file)
    if not match:
        raise ValueError(f"Unexpected filename {df_file}, did not match regexp")
    else:
        n = int(match.group(1))
        mriot = match.group(2)
        df = process_impact_df(df_file, n, mriot)
        dfs_impact.append(df.reset_index())

meta_df_imp = pd.concat(dfs_impact, ignore_index=True)
dfs_res = []


for res_file in snakemake.input.results:
    match = regex_results.search(res_file)
    if not match:
        raise ValueError(f"Unexpected filename {res_file}, did not match regexp")
    else:
        n = int(match.group(1))
        mriot = match.group(2)
        simtype = match.group(3)
        var = match.group(4)
        df = process_res_df(res_file, n, mriot, simtype, var)
        dfs_res.append(df.reset_index())

meta_df_res = pd.concat(dfs_res, ignore_index=True)

meta_df_imp.set_index(["sample","MRIOT","region","sector"],inplace=True)
meta_df_res.set_index(["sample","MRIOT","region","sector"],inplace=True)

meta_df = meta_df_imp.join(meta_df_res)
meta_df.to_parquet(pathlib.Path(snakemake.output.meta_results))
