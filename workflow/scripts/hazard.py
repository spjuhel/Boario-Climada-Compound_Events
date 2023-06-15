import sys
from climada.util.api_client import Client
from climada.engine.impact_calc import ImpactCalc
from climada.entity import ImpfSetTropCyclone
from climada.hazard import Hazard
import numpy as np
from datetime import datetime
import logging
import pickle as pkl

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

logger.info("Starting")
hazfilename_list = snakemake.params.hazfiles
year_start = snakemake.params.year_start
year_end = snakemake.params.year_end

client = Client()
# getting the right id for he impact function in our exposures
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
