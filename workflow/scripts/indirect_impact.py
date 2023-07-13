import logging
import pathlib
import pandas as pd
import pickle as pkl
from boario.extended_models import ARIOPsiModel
from boario.event import EventKapitalRebuild
from boario.simulation import Simulation

logFormatter = logging.Formatter(
    "%(asctime)s [%(levelname)-5.5s] %(name)s %(message)s", datefmt="%H:%M:%S"
)
scriptLogger = logging.getLogger("Supchain from Impact")
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
scriptLogger.addHandler(consoleHandler)
scriptLogger.setLevel(logging.INFO)


def simulate(
    *,
    mriot,
    inv_dict,
    k_vector,
    rebuild_tau,
    order_type,
    alpha_tau,
    monetary_factor,
    psi_param,
    inventory_restoration_tau,
    df_impact,
    reb_sect,
    duration,
    output,
    vars_to_save,
    event_id=None,
):
    model = ARIOPsiModel(
        mriot,
        order_type=order_type,
        alpha_tau=alpha_tau,
        rebuild_tau=rebuild_tau,
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
        n_temporal_units_to_sim=int(df_impact.index.max()) + 365,
        separate_sims=False,
    )

    if event_id:
        impact_row = df_impact.iloc[event_id]
        events_list = [
            EventKapitalRebuild.from_series(
                impact=impact_row,
                rebuilding_sectors=reb_sect,
                rebuild_tau=rebuild_tau,
                rebuilding_factor=1.0,
                households_impact=[],
                occurrence=impact_row.name+1,
                duration=duration,
                # TODO: it is assumed that exposure/impacts are expressed with full
                # values, deal with cases where CLIMADA exposure is expressed in K,
                # M or Bn
                event_monetary_factor=1.0,
            )
        ]
    else:
        events_list = [
            EventKapitalRebuild.from_series(
                impact=ev[1],
                rebuilding_sectors=reb_sect,
                rebuild_tau=rebuild_tau,
                rebuilding_factor=1.0,
                households_impact=[],
                occurrence=ev[0]+1,
                duration=duration,
                # TODO: it is assumed that exposure/impacts are expressed with full
                # values, deal with cases where CLIMADA exposure is expressed in K,
                # M or Bn
                event_monetary_factor=1.0,
            )
            for ev in df_impact.iterrows()
        ]

    sim.add_events(events_list)
    sim.loop(progress=False)
    for var_name in vars_to_save:
        var_value = getattr(sim, var_name)
        var_value.to_parquet(pathlib.Path(output) / f"{var_name}.parquet")


with open(snakemake.input.supchain, "rb") as f:
    supchain = pkl.load(f)

with open(snakemake.input.mriot_file, "rb") as f:
    mriot = pkl.load(f)

sectors_df = pd.read_parquet(snakemake.input.sectors_df)
reb_sect = sectors_df.loc[
    sectors_df.rebuilding_factor > 0, "rebuilding_factor"
].to_dict()
inv_dict = sectors_df.loc[:, "inventory_size"].to_dict()

event_id = snakemake.wildcards.event_id
if event_id == "aggregated":
    event_id = None
else:
    event_id = int(event_id)

df_impact=pd.read_parquet(snakemake.input.df_impact)

k_vector = (supchain.secs_exp / supchain.conversion_factor())

simulate(
    rebuild_tau=snakemake.params.rebuild_tau,
    order_type=snakemake.params.order_type,
    alpha_tau=snakemake.params.alpha_tau,
    monetary_factor=snakemake.params.monetary_factor,
    psi_param=snakemake.params.psi,
    inventory_restoration_tau=snakemake.params.inv_tau,
    duration=snakemake.params.duration,
    mriot=mriot,
    inv_dict=inv_dict,
    k_vector=k_vector,
    df_impact=df_impact,
    reb_sect=reb_sect,
    output=snakemake.output.path,
    event_id=event_id,
    vars_to_save=snakemake.config["variables"]
)
