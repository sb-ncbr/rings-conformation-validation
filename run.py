import time
import traceback
import logging
from datetime import timedelta
from workflow.analysis.conformer_assignment import assign_conformations
from workflow.analysis.coverage import analyse_coverage
from workflow.data_acquisition.rings_extraction import extract_rings
from utils.config import load_config
from utils.logging import setup_logging
from workflow.data_acquisition.download import download_ccp4_with_retries
from workflow.postprocessing.build_web_dataset import create_data_for_web
from workflow.preprocessing.rings_filtering import filter_rings
from workflow.utils.timing import ExecutionTimer, timed_step


def main():
    start = time.perf_counter()

    logger = logging.getLogger(__name__)
    config = load_config()
    setup_logging(config.log_dir)

    logger.info("Workflow started")

    try:
        timer = ExecutionTimer()

        config.main_dir.mkdir(parents=True, exist_ok=True)

        with timed_step(timer, "Extract rings"):
            extract_rings(config.pdb_dir, config.main_dir)

        with timed_step(timer, "Filter rings"):
            filter_rings(config.main_dir, config.ccd, config.patterns_dir, config.state_dir)

        with timed_step(timer, "Assign conformations"):
            assign_conformations(config.main_dir)
        
        with timed_step(timer, "Download density maps"):
            download_ccp4_with_retries(config.density_dir, config.main_dir)

        with timed_step(timer, "Electron density coverage"):
            analyse_coverage(config.main_dir, config.density_dir, False, False)

        with timed_step(timer, "Build web dataset"):
            create_data_for_web(config.output_dir, config.main_dir, config.pdb_dir, config.methods_info)

        logger.info("Workflow completed successfully")

    except Exception:
        logger.error("Workflow failed")
        logger.error(traceback.format_exc())
        raise

    finally:
        logger.info("Execution time summary:")
        timer.report()

        elapsed = time.perf_counter() - start
        formatted = str(timedelta(seconds=elapsed))

        logger.info(
            "Total execution time: %s",
            formatted
        )



if __name__ == "__main__":
    main()