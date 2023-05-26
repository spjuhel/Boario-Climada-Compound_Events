import pickle as pkl

import pandas as pd

from boario.extended_models import ARIOPsiModel
from boario.event import *
from boario.simulation import Simulation


def simulate(
    separate,
    event_number,
    output_conso,
    output_prod,
    rebuild_tau,
    duration,
    inv_dict,
    reb_sect,
    k_vector,
    mriot,
    order_type="alt",
    alpha_tau=365,
    main_inv_dur=90,
    monetary_factor=1000000,
    psi_param=0.80,
    inventory_restoration_tau=60,
):
    model = ARIOPsiModel(
        mriot,
        order_type=order_type,
        alpha_tau=alpha_tau,
        main_inv_dur=main_inv_dur,
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
        n_temporal_units_to_sim=df_impact.index.max() + 365,
        separate_sims=False,
    )

    if separate:
        events_list = [
            EventKapitalRebuild.from_series(
                impact=ev[1],
                rebuilding_sectors=reb_sect,
                rebuild_tau=rebuild_tau,
                rebuilding_factor=1.0,
                households_impact=[],
                occurrence=ev[0],  # check this
                duration=duration,
                event_monetary_factor=1.0,
            )
            for ev in df_impact.iloc[[event_number]].iterrows()
        ]

    else:
        events_list = [
            EventKapitalRebuild.from_series(
                impact=ev[1],
                rebuilding_sectors=reb_sect,
                rebuild_tau=rebuild_tau,
                rebuilding_factor=1.0,
                households_impact=[],
                occurrence=ev[0],  # check this
                duration=duration,
                event_monetary_factor=1.0,
            )
            for ev in df_impact.iterrows()
        ]
        sim.add_events(events_list)
        sim.loop()
        aggregated_prod_sim = sim.production_realised.copy()
        aggregated_prod_sim.to_parquet(output_prod)

    aggregated_conso_sim = sim.final_demand_not_met.copy()
    aggregated_conso_sim.to_parquet(output_conso)


with open(snakemake.input.mriot, "rb") as f:
    mriot_arg = pkl.load(f)

k_vector_arg = pd.read_parquet(snakemake.input.k_vector)
df_impact = pd.read_parquet(snakemake.input.df_impact)
sectors_df = pd.read_parquet(snakemake.input.sectors_df)

reb_sect_arg = sectors_df.loc[
    sectors_df.rebuilding_factor > 0, "rebuilding_factor"
].to_dict()
inv_dict_arg = sectors_df.loc[:, "inventory_size"].to_dict()

if snakemake.params.separate:
    simulate(
        separate=snakemake.params.separate,
        event_number=snakemake.wildcards.i,
        output_conso=snakemake.output.conso,
        output_prod=snakemake.output.prod,
        rebuild_tau=snakemake.params.rebuild_tau,
        duration=snakemake.params.duration,
        inv_dict=inv_dict_arg,
        reb_sect=reb_sect_arg,
        k_vector=k_vector_arg,
        mriot=mriot_arg,
        order_type=snakemake.params.order_type,
        alpha_tau=snakemake.params.alpha_tau,
        main_inv_dur=snakemake.params.main_inv_dur,
        monetary_factor=1000000,
        psi_param=snakemake.params.psi,
        inventory_restoration_tau=snakemake.params.inv_tau,
    )
else:
    simulate(
        separate=snakemake.params.separate,
        event_number=0,
        output_conso=snakemake.output.conso,
        output_prod=snakemake.output.prod,
        rebuild_tau=snakemake.params.rebuild_tau,
        duration=snakemake.params.duration,
        inv_dict=inv_dict_arg,
        reb_sect=reb_sect_arg,
        k_vector=k_vector_arg,
        mriot=mriot_arg,
        order_type=snakemake.params.order_type,
        alpha_tau=snakemake.params.alpha_tau,
        main_inv_dur=snakemake.params.main_inv_dur,
        monetary_factor=1000000,
        psi_param=snakemake.params.psi,
        inventory_restoration_tau=snakemake.params.inv_tau,
    )
