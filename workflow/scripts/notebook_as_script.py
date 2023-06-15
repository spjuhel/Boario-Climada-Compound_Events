from climada.util.api_client import Client
from climada_petals.engine import SupplyChain
from climada.entity import ImpfSetTropCyclone, ImpactFuncSet
from climada.engine.impact_calc import ImpactCalc
from climada.hazard import Hazard
import numpy as np
import pathlib
from functools import reduce
from datetime import datetime

def get_impf_id(cnt):
    for basin, iso_list in ImpfSetTropCyclone.get_countries_per_region()[2].items():
        for iso in iso_list:
            if iso == cnt:
                return basin, ImpfSetTropCyclone.get_countries_per_region()[1][basin]

if __name__ == "main":

    client = Client()
    # getting the right id for he impact funcion in our exposures
    tc_wp_1 = client.get_hazard('tropical_cyclone', name='tropical_cyclone_0synth_tracks_150arcsec_genesis_WP_1980_2020')
    tc_na = client.get_hazard('tropical_cyclone', name='tropical_cyclone_0synth_tracks_150arcsec_genesis_NA_1980_2020')

    tc_wp = Hazard.concat([tc_wp_1,tc_na])

    assets = client.get_litpop()


    impf_TC = ImpfSetTropCyclone.from_calibrated_regional_ImpfSet()

    assets.gdf['impf_TC'] = 1

    for cnt in np.unique(assets.gdf.region_id):
        assets.gdf.loc[assets.gdf['region_id']==cnt, 'impf_TC'] = get_impf_id(int(cnt))[1]

    # impact calculation, already selecting the date
    imp_calc = ImpactCalc(assets, impf_TC, tc_wp)

    # Create a datetime object
    dt_start = datetime(2012, 1, 1)
    dt_end = datetime(2017, 1, 1)

    # Get the proleptic Gregorian ordinal value
    ordinal_value_start = dt_start.toordinal()
    ordinal_value_end = dt_end.toordinal()

    direct_impact = imp_calc.impact().select(dates=(ordinal_value_start,ordinal_value_end))

    supchain = SupplyChain.from_mriot(mriot_type='WIOD16', mriot_year=2011)

    impacted_secs = supchain.mriot.get_sectors().tolist()

    supchain.calc_secs_exp_imp_shock(assets, direct_impact, impacted_secs)
    supchain.calc_direct_production_impacts()

    dir_prod_impt_mat = supchain.dir_prod_impt_mat

    supchain_impacted = supchain.dir_prod_impt_mat.loc[:, (supchain.dir_prod_impt_mat > 0).any()]
    supchain_impacted = supchain_impacted.loc[(supchain_impacted > 0).any(axis=1), :]

    matching_indices = np.where(np.isin(tc_wp.event_id, supchain_impacted.index))[0]

    event_date = tc.date[matching_indices]
    supchain_impacted["date"] = event_date

    supchain_impacted.fillna(0, inplace=True)

    save_path.mkdir(exist_ok=True)

    supchain.calc_indirect_production_impacts("boario_aggregated")
    supchain.indir_prod_impt_mat.to_parquet(save_path / "aggregated.parquet")

    disaggregated_app_res = supchain.calc_indirect_production_impacts("boario_separated")
    disaggregated_app_res_red = reduce(lambda x, y: x.add(y-y.loc[0], fill_value=0), disaggregated_app_res)
    disaggregated_app_res_red.to_parquet(save_path / "separated.parquet")
