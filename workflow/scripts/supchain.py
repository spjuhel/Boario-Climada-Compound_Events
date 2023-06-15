from climada_petals.engine import SupplyChain
import pandas as pd
import logging
import pickle as pkl

logFormatter = logging.Formatter(
    "%(asctime)s [%(levelname)-5.5s] %(name)s %(message)s", datefmt="%H:%M:%S"
)
scriptLogger = logging.getLogger("Supchain from Impact")
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
scriptLogger.addHandler(consoleHandler)
scriptLogger.setLevel(logging.INFO)

def load_sectors_aggreg(mrio_name,sectors_common_aggreg):
    mrio_name = mrio_name.casefold()
    if "eora" in mrio_name:
        return sectors_common_aggreg["eora26_without_reexport_to_common_aggreg"]
    elif "euregio" in mrio_name:
        return sectors_common_aggreg["euregio_to_common_aggreg"]
    elif "exio" in mrio_name:
        return sectors_common_aggreg["exiobase_full_to_common_aggreg"]
    elif "oecd" in mrio_name:
        return sectors_common_aggreg["icio2021_reworked_to_common_aggreg"]
    else:
        raise ValueError(f"Invalid MRIO name: {mrio_name}")


with open(snakemake.input.impact_file,"rb") as f:
    direct_impact = pkl.load(f)

with open(snakemake.input.mriot_file,"rb") as f:
    mriot = pkl.load(f)

with open(snakemake.input.assets,"rb") as f:
    assets = pkl.load(f)

sectors_common_aggreg = {
    sheet_name: pd.read_excel(snakemake.params.sector_common_aggreg, sheet_name=sheet_name, index_col=0)
for sheet_name in [
        "eora26_without_reexport_to_common_aggreg",
        "euregio_to_common_aggreg",
        "exiobase_full_to_common_aggreg",
        "icio2021_reworked_to_common_aggreg",
]
}

df_aggreg = load_sectors_aggreg(mriot.name, sectors_common_aggreg)

mriot.rename_sectors(df_aggreg["new sector"].to_dict())
mriot.aggregate_duplicates()

mriot.meta.change_meta("name", snakemake.wildcards.mriot_name)

if "EXIOBASE3" in mriot.meta.name:
    reg = mriot.get_regions()
    agg_regions = [rg if rg not in ['WA', 'WE',  'WF',  'WL',  'WM'] else "ROW" for rg in reg]
    mriot = mriot.aggregate(region_agg = agg_regions)

mriot.reset_all_full()
mriot.calc_all()

supchain = SupplyChain(mriot)
impacted_secs = supchain.mriot.get_sectors().tolist()

# First pass
supchain.calc_shock_to_sectors(assets, direct_impact, impacted_secs)

supchain_impacted = supchain.secs_imp.loc[:, (supchain.secs_imp > 0).any()]
supchain_impacted = supchain_impacted.loc[(supchain_impacted > snakemake.params.direct_dmg_floor).any(axis=1), :]

# Second pass to eliminate events with dmg < snakemake.params.direct_dmg_floor
direct_impact = direct_impact.select(event_ids=supchain_impacted.index)
supchain.calc_shock_to_sectors(assets, direct_impact, impacted_secs)

df_impact = pd.DataFrame(
        supchain.secs_imp.values,
        columns=supchain.secs_imp.columns,
    index = supchain.events_date
    )
df_impact.index = df_impact.index - df_impact.index.min()

with open(snakemake.output.mriot_save,"wb") as f:
    pkl.dump(mriot,f)

with open(snakemake.output.save_path,"wb") as f:
    pkl.dump(supchain,f)

df_impact.to_parquet(snakemake.output.df_impact)
