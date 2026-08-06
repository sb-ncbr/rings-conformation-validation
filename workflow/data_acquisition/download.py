from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import logging
from typing import Set
import requests
import gzip
from pathlib import Path
from workflow.utils.helpers import get_pdb_names_filtered


logger = logging.getLogger(__name__)

CCP4_NOT_FOUND_CACHE = Path("ccp4_not_found.json")
CCP4_FAILED_CACHE = Path("ccp4_failed.json")


def load_json(filename):
    if filename.exists():
        with open(filename) as f:
            return set(json.load(f))
    return set()


def save_to_json(data: Set, filename):
    with open(filename, "w") as f:
        json.dump(sorted(data), f, indent=2)


def download_and_gzip(pdb_id, out_dir: Path, timeout=20, chunk_size=8192):
    name = pdb_id.lower()
    BASE_URL = "https://www.ebi.ac.uk/pdbe/coordinates/files/"
    url = "".join([BASE_URL, name, ".ccp4"])
    out_path = out_dir / f"{name}.ccp4.gz"

    try:
        with requests.get(url, stream=True, timeout=timeout) as r:
            if r.status_code != 200:
                
                if r.status_code == 404:
                    # no data available for that structure -> ok
                    return name, True, "not_found"
                
                logger.error(
                    f"[{name}] "
                    f"(status={r.status_code}, reason={r.reason}) | {url}"
                )
                
                return name, False, f"HTTP {r.status_code}"

            with gzip.open(out_path, "wb") as f_out:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    if chunk:  # filter keep-alive chunks
                        f_out.write(chunk)
        logger.info(f"Successfully downloaded and zipped: {out_path.name}")
        return name, True, None

    except Exception as e:
        return name, False, str(e)


def download_many(names: Set[str], not_found: Set, out_dir: Path, max_workers: int):
    out_dir.mkdir(parents=True, exist_ok=True)

    results = {
        "success": 0,
        "not_found": 0,
        "failed": 0
    }

    # not_found = set()
    failed = set()

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(download_and_gzip, n, out_dir) for n in names]

        for future in as_completed(futures):
            name, ok, err = future.result()

            if ok: 
                if err == "not_found":
                    not_found.add(name.lower())
                    results["not_found"] += 1
                    
                else:
                    results["success"] += 1
            else:
                failed.add(name.lower())
                results["failed"] += 1

    save_to_json(not_found, CCP4_NOT_FOUND_CACHE)
    save_to_json(failed, CCP4_FAILED_CACHE)

    return results


def download_density_maps(ccp4_dir: Path, main_dir: Path):
    logger.info("Downloading density maps...")

    names = get_pdb_names_filtered(main_dir)
    ccp4_dir.mkdir(parents=True, exist_ok=True)

    not_found = load_json(CCP4_NOT_FOUND_CACHE)
    prev_failed = load_json(CCP4_FAILED_CACHE)
    logger.info(f'Previously failed: {len(prev_failed)}')
    logger.info(f'Not found: {len(not_found)}')

    names = [
        pdb
        for pdb in names
        if pdb.lower() not in not_found
        and not (
        (path := ccp4_dir / f"{pdb.lower()}.ccp4.gz").exists()
        and path.stat().st_size > 0
    )
    ]

    logger.info(f"{len(names)} files will be processed...")


    results = download_many(names, not_found, out_dir=ccp4_dir, max_workers=32)

    logger.info(f"Downloaded: {results["success"]}")
    logger.info(f"Not found: {results["not_found"]}")
    logger.info(f"Failed now: {results["failed"]}")

    return results["failed"]


def download_ccp4_with_retries(density_dir: Path, main_dir: Path):
    max_retries = 5

    for attempt in range(max_retries):
        failed = download_density_maps(density_dir, main_dir)

        if failed == 0:
            break

        logger.warning(f"Download attempt {attempt + 1}/{max_retries}: {failed} failed")

    else:
        raise RuntimeError(f"Failed to download {failed} density maps after {max_retries} attempts")
    

