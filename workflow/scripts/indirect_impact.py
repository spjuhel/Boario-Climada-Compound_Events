import sys
import logging
import pathlib
import pandas as pd
import pickle as pkl
from boario.extended_models import ARIOPsiModel
from boario.event import EventKapitalRebuild
from boario.simulation import Simulation
from boario import logger as lg

import traceback

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

def load_sectors_aggreg(mrio_name,sectors_common_aggreg):
    mrio_name = mrio_name.casefold()
    if "eora" in mrio_name:
        return sectors_common_aggreg["eora26_without_reexport_to_common_aggreg"]
    elif "euregio" in mrio_name:
        return sectors_common_aggreg["euregio_to_common_aggreg"]
    elif "exio" in mrio_name:
        return sectors_common_aggreg["exiobase_full_to_common_aggreg"]
    elif "oecd" in mrio_name:
        return sectors_common_aggreg["icio2021_to_common_aggreg"]
    else:
        raise ValueError(f"Invalid MRIO name: {mrio_name}")

sys.excepthook = handle_exception

logFormatter = logging.Formatter(
    "%(asctime)s [%(levelname)-5.5s] %(name)s %(message)s", datefmt="%H:%M:%S"
)
ch.setFormatter(logFormatter)
fh.setFormatter(logFormatter)
logger.addHandler(ch)
logger.addHandler(fh)
lg.addHandler(fh)
logger.setLevel(logging.INFO)
lg.setLevel(logging.INFO)

def simulate(
    *,
    mriot,
    inv_dict,
    k_vector,
    #rebuild_tau,
    order_type,
    alpha_tau,
    alpha_max,
    monetary_factor,
    psi_param,
    inventory_restoration_tau,
    df_impact,
    meta_df_impact,
    reb_sect,
    duration,
    output,
    vars_to_save,
    step_to_eq,
    event_id=None,
):
    model = ARIOPsiModel(
        mriot,
        order_type=order_type,
        alpha_tau=alpha_tau,
        alpha_max=alpha_max,
        monetary_factor=monetary_factor,
        inventory_dict=inv_dict,
        productive_capital_vector=k_vector,
        psi_param=psi_param,
        inventory_restoration_tau=inventory_restoration_tau,
    )

    # run simulation up to one year after the last event
    sim = Simulation(
        model,
        register_stocks=False,
        n_temporal_units_to_sim=int(df_impact.index.max()) + step_to_eq,
        separate_sims=False,
    )

    if event_id is not None:
        impact_row = df_impact.iloc[event_id]
        meta_impact_row = meta_df_impact.iloc[event_id]
        events_list = [
            EventKapitalRebuild.from_series(
                impact=impact_row,
                rebuilding_sectors=reb_sect,
                rebuild_tau=meta_impact_row["recovery_duration"],
                rebuilding_factor=1.0,
                households_impact=[],
                occurrence=impact_row.name+1,
                duration=duration,
                event_monetary_factor=1.0,
            )
        ]
    else:
        events_list = [
            EventKapitalRebuild.from_series(
                impact=ev[1],
                rebuilding_sectors=reb_sect,
                rebuild_tau=meta_df_impact.loc[ev[0],"recovery_duration"],
                rebuilding_factor=1.0,
                households_impact=[],
                occurrence=ev[0]+1,
                duration=duration,
                event_monetary_factor=1.0,
            )
            for ev in df_impact.iterrows()
        ]

    sim.add_events(events_list)
    sim.loop(progress=False)
    for var_name in vars_to_save:
        var_value = getattr(sim, var_name)
        var_value.to_parquet(pathlib.Path(output) / f"{var_name}.parquet")

logger.info("Starting")

with open(snakemake.input.supchain, "rb") as f:
    supchain = pkl.load(f)

with open(snakemake.input.mriot_file, "rb") as f:
    mriot = pkl.load(f)


if snakemake.config["aggregate_mriot"]:
    sectors_common_aggreg = {
        sheet_name: pd.read_excel(snakemake.params.sector_common_aggreg, sheet_name=sheet_name, index_col=0)
        for sheet_name in [
                "eora26_without_reexport_to_common_aggreg",
                "euregio_to_common_aggreg",
                "exiobase_full_to_common_aggreg",
                "icio2021_to_common_aggreg",
        ]
    }
    df_aggreg = load_sectors_aggreg(mriot.name, sectors_common_aggreg)

    mriot.rename_sectors(df_aggreg["new sector"].to_dict())
    mriot.aggregate_duplicates()

sectors_df = pd.read_parquet(snakemake.input.sectors_df)
reb_sect = sectors_df.loc[
    sectors_df.rebuilding_factor > 0, "rebuilding_factor"
].to_dict()
inv_dict = sectors_df.loc[:, "inventory_size"].to_dict()
inv_tau = sectors_df.loc[:, "inventory_tau"].to_dict()

event_id = snakemake.wildcards.event_id
if event_id == "aggregated":
    event_id = None
else:
    event_id = int(event_id)

meta_df_impact=pd.read_parquet(snakemake.input.meta_df_impact)

df_impact=pd.read_parquet(snakemake.input.df_impact)
k_vector = (supchain.secs_exp / supchain.conversion_factor())

if snakemake.config["aggregate_mriot"]:
    df_impact = df_impact.rename(df_aggreg["new sector"].to_dict(), axis=1, level=1)
    df_impact = df_impact.T.groupby(['region','sector']).sum().T

    k_vector = k_vector.rename(df_aggreg["new sector"].to_dict(), axis=1, level=1)
    k_vector = k_vector.T.groupby(['region','sector']).sum().T

simulate(
    #rebuild_tau=snakemake.params.rebuild_tau,
    order_type=snakemake.params.order_type,
    alpha_tau=snakemake.params.alpha_tau,
    alpha_max=snakemake.params.alpha_max,
    monetary_factor=snakemake.params.monetary_factor,
    psi_param=snakemake.params.psi,
    inventory_restoration_tau=inv_tau,
    duration=snakemake.params.duration,
    mriot=mriot,
    inv_dict=inv_dict,
    k_vector=k_vector,
    df_impact=df_impact,
    meta_df_impact=meta_df_impact,
    reb_sect=reb_sect,
    output=snakemake.output.path,
    step_to_eq=snakemake.params.step_to_eq,
    event_id=event_id,
    vars_to_save=snakemake.params.variables
)

logger.info("Finished")
