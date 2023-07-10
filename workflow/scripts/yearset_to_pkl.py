import sys
import logging
import pickle as pkl

from climada.engine import Impact
from scipy.sparse import csr_matrix

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

def read_csv_npz(input_csv_path, sparse_csr_path):
    impact = Impact.from_csv(input_csv_path)
    # Read the sparse CSR matrix from a file
    impact.imp_mat = Impact.read_sparse_csr(sparse_csr_path)
    return impact

logger.info("Starting")

direct_impact = read_csv_npz(snakemake.input.csv, snakemake.input.npz)

with open(snakemake.output.pkl,"wb") as f:
    pkl.dump(direct_impact,f)
logger.info("Finished")
