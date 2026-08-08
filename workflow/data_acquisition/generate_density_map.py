from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
from pathlib import Path
import time
import pandas as pd
import subprocess
import gzip
import shutil
import os

from workflow.models.Ring import Ring

# from workflow.utils.helpers import get_pdb_names_filtered

logger = logging.getLogger(__name__)

def get_pdb_names_filtered(main_dir: Path):
    csv_files = []
    for ring in Ring:
        csv_files.append(os.path.join(str(main_dir), ring.name.lower(), f"filtered_patterns_{ring.name.lower()}.csv"))

    unique_ids = set()

    for f in csv_files:
        for chunk in pd.read_csv(f, usecols=["ParentId"], chunksize=100_000):
            # unique_ids.update(chunk["ParentId"].dropna())
            unique_ids.update(chunk["ParentId"].dropna().astype(str).str[-12:]) # changed to last 12 chars
            # 0k_pdb_000010kt_structures__pdb_000010kt

    logger.info(f"Unique pdb ids across all filtered rings: {len(unique_ids)}")
    return unique_ids



def make_ccp4_gz(args):
    cif_path, out_path = args

    cif_path = Path(cif_path)
    out_path = Path(out_path)
    output_gz = Path(str(out_path) + ".gz")

    pdb_id = out_path.stem

    start_total = time.perf_counter()

    try:
        
        start = time.perf_counter()
        # logger.info(f"Generating CCP4: {cif_path.name}")

        subprocess.run(
            ["gemmi", "sf2map", str(cif_path), str(out_path)],
            check=True,
            capture_output=True,
            text=True
        )
        map_time = time.perf_counter() - start

        
        # logger.info(f"Compressing: {out_path.name}")
        start = time.perf_counter()

        with open(out_path, "rb") as f_in:
            with gzip.open(output_gz, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        compress_time = time.perf_counter() - start

        out_path.unlink()

        total_time = time.perf_counter() - start_total
        # logger.info(f"Finished: {output_gz.name}")

        # logger.info(
        #     f"{pdb_id} | "
        #     f"map generated in {map_time:.2f}s | "
        #     f"compressed in {compress_time:.2f}s | "
        #     f"total {total_time:.2f}s"
        # )


        return True, {  
            "pdb_id": pdb_id,
            "map_time": map_time,
            "gzip_time": compress_time,
            "total_time": total_time
            }
    
    except subprocess.CalledProcessError as e:
        logger.error(
            f"Gemmi failed for {pdb_id}: "
            f"{e.stderr.strip() if e.stderr else str(e)}"
        )

        # cleanup partial file if created
        if out_path.exists():
            out_path.unlink()

        return False, None

    except Exception as e:
        logger.exception(f"Failed processing {pdb_id}: {e}")

        # cleanup partial files
        for p in (out_path, output_gz):
            if p.exists():
                p.unlink()

        return False, None


def run_ccp4_generation(main_dir: Path, pdb_dir: Path,  density_dir: Path):
    density_dir.mkdir(parents=True, exist_ok=True)
    pdbs = get_pdb_names_filtered(main_dir)
    jobs = []
    for pdb_id in pdbs:
        input_path = pdb_dir / pdb_id[9:11] / pdb_id / "validation_reports" / f"{pdb_id}_validation_2fo-fc_map_coef.cif.gz"
        output_path = density_dir / f"{pdb_id}.ccp4"
        gz_output = Path(str(output_path) + '.gz')
        if gz_output.exists() and gz_output.stat().st_size > 0:
            continue
        if input_path.exists():
            jobs.append((input_path, output_path))

    n_workers = os.cpu_count()

    results = []

    total = len(jobs)
    completed = 0

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = [executor.submit(make_ccp4_gz, job) for job in jobs]

        for future in as_completed(futures):
            result = future.result()
            completed += 1

            success, info = result

            if success:
                logger.info(
                    f"[{completed}/{total}] {info['pdb_id']} | "
                    f"map {info['map_time']:.2f}s | "
                    f"gzip {info['gzip_time']:.2f}s | "
                    f"total {info['total_time']:.2f}s"
                )

    # with ProcessPoolExecutor(max_workers=n_workers) as executor:
    #     for result in executor.map(make_ccp4_gz, jobs):
    #         results.append(result)

    # failed = [r for r in results if not r[1]]

    # logger.info(f"Successfully finished: {len(results)-len(failed)}/{len(results)}")
    # logger.info(f"Failed: {len(failed)}")