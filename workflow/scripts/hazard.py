import sys
import numpy as np
from datetime import datetime
import logging
import pickle as pkl
import pycountry

from climada.util.api_client import Client
from climada.engine.impact_calc import ImpactCalc
from climada.entity import ImpfSetTropCyclone
from climada.hazard import Hazard
import climada.util.dates_times as u_dt
from climada.util import yearsets, coordinates
from climada.util import yearsets

logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

logger = logging.getLogger(__name__)

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


# Install exception handler
sys.excepthook = handle_exception

logFormatter = logging.Formatter(
    "%(asctime)s [%(levelname)-5.5s] %(name)s %(message)s", datefmt="%H:%M:%S"
)
scriptLogger = logger
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
scriptLogger.addHandler(consoleHandler)
scriptLogger.setLevel(logging.INFO)

def get_impf_id(cnt):
    for basin, iso_list in ImpfSetTropCyclone.get_countries_per_region()[2].items():
        for iso in iso_list:
            if iso == cnt:
                return basin, ImpfSetTropCyclone.get_countries_per_region()[1][basin]

def get_hazard_from_api(hazfilename_list:"list[(str,str)]", client):
    return [client.get_hazard(haztype, name=hazname) for haztype, hazname in hazfilename_list]

def year_date_event_in_sample(years, dates, sampling_vec):
    """
    Change the year for the sampled events

    Parameters
    ----------
    years : list[int]
        Years of sampled events (length equal to length of sampling_vec)
    dates : list[dates]
        List of dates in ordinal format for the whole event set
    sampling_vec : list[np.array]
        Array of ids (row index) of selected events per year.

    Raises
    ------
    ValueError
        length of years must equal length of sampling_vec

    Returns
    -------
    list[dates]
        List with dates in ordinal format of the sampled events.

    """
    if len(years) != len(sampling_vec):
        raise ValueError("The number of years is different from the length" +
                         "of the sampling vector")
    def change_year(old_date, year):
        old_date = u_dt.date_to_str(old_date)
        if old_date[5:7] == '02' and old_date[8:10] == '29':
            old_date = old_date.replace('29', '28')
        new_date = old_date.replace(old_date[0:4], str(year))
        return u_dt.str_to_date(new_date)

    return [
        change_year(date, year)
        for year, events in zip(years, sampling_vec)
        for date in np.array(dates)[events]
        ]

def make_yearset(impact_dict, n_years):
    """
    Parameters
    ----------
    impact_dict : dict of Impact
        dictionary of impacts
    n_years : number of years

    Returns
    -------
    imp : dict

    """
    years = np.arange(1, n_years+1)
    yearset_dict = {}
    for impact in impact_dict: # we do this per basin as the frequencies are different
        lam = np.sum(impact_dict[impact].frequency)
        lam = np.round(lam, 10)
        events_per_year = yearsets.sample_from_poisson(len(years), lam) # number of events to pick for each year
        sampling_vect = yearsets.sample_events(events_per_year, impact_dict[impact].frequency) # vector containing the the event id per year
        yearset = copy.deepcopy(impact_dict[impact]) # new impact object
        yearset.date = year_date_event_in_sample(years,impact_dict[impact].date,sampling_vect)
        sampling_vect = np.concatenate(sampling_vect)
        yearset.event_name = list(np.array(impact_dict[impact].event_name)[sampling_vect])
        yearset.at_event = impact_dict[impact].at_event[sampling_vect]
        yearset.event_id = np.arange(1, len(yearset.event_name)+1)
        yearset.frequency = np.ones(len(yearset.event_id))/n_years

        yearset.imp_mat = impact_dict[impact].imp_mat[sampling_vect, :]
        yearset_dict[impact] = yearset
    return yearset_dict

# create n years of impacts by sampling events based on the expected number of events. This is
# done per basin
def yearset_impact(direct_impact_wp, direct_impact_na, countries_num, n_years):

    impact_dict = {'wp':direct_impact_wp, 'na':direct_impact_na}
    ys = make_yearset(impact_dict, n_years) # we get n_years for each basin separatly
    imp = Impact()
    imp.frequency = np.hstack([ys['wp'].frequency, ys['na'].frequency])
    imp.date = np.hstack([ys['wp'].date, ys['na'].date])
    imp.event_name = np.hstack([ys['wp'].event_name, ys['na'].event_name])

    imp.event_id = np.hstack([ys['wp'].event_id, (ys['na'].event_id+np.max(ys['wp'].event_id))])
    imp.coord_exp = direct_impact_wp.coord_exp
    imp.imp_mat = sparse.vstack([ys['wp'].imp_mat,ys['na'].imp_mat])


    imp.at_event, imp.eai_exp, imp.aai_agg = ImpactCalc.risk_metrics(imp.imp_mat, imp.frequency)
    imp = imp.select(event_ids = imp.event_id[imp.at_event>1e8])
    imp.tag = impact_dict['wp'].tag
    country_matrices = {country: np.array([imp.imp_mat[:, countries_num == countries[country]].sum(axis=1).flatten()])
                     for country in countries}

    years = [datetime.date.fromordinal(date).year for date in imp.date]
    country_year_impact = {year: {country: country_matrices[country][0][0][years==year].sum()\
                              for country in countries} for year in np.unique(years)}

    n_countries_per_year = {year:np.sum([np.array(list(country_year_impact[year].values()))>1e8]) for year in years}


    return (imp, country_year_impact, n_countries_per_year, imp.aai_agg)

logger.info("Starting")
hazfilename_list = snakemake.params.hazfiles
year_start = snakemake.params.year_start
year_end = snakemake.params.year_end

client = Client()
# getting the right id for the impact function in our exposures
hazards = get_hazard_from_api(hazfilename_list=hazfilename_list, client=client)

hazards_concat = Hazard.concat(hazards)
assets = client.get_litpop()

with open(snakemake.output.assets,"wb") as f:
    pkl.dump(assets,f)

# TODO : make this less arbitratry/hardcoded
impf_TC = ImpfSetTropCyclone.from_calibrated_regional_ImpfSet()

assets.gdf['impf_TC'] = 1

for cnt in np.unique(assets.gdf.region_id):
    assets.gdf.loc[assets.gdf['region_id']==cnt, 'impf_TC'] = get_impf_id(int(cnt))[1]

# impact calculation, already selecting the date
imp_calc = ImpactCalc(assets, impf_TC, hazards_concat)

# Create a datetime object
dt_start = datetime(year_start, 1, 1)
dt_end = datetime(year_end, 1, 1)

# Get the proleptic Gregorian ordinal value
ordinal_value_start = dt_start.toordinal()
ordinal_value_end = dt_end.toordinal()

direct_impact = imp_calc.impact().select(dates=(ordinal_value_start,ordinal_value_end))
with open(snakemake.output.save_path,"wb") as f:
    pkl.dump(direct_impact,f)
logger.info("Finished")
