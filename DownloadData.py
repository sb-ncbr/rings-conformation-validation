#!/usr/bin/env python3

"""
Python script by which you can download Onedata space or directory with all its content. Downloaded directories and files have to be shared in Onedata.
"""

import argparse
import os
import sys
import random
import re
import threading
import queue
from multiprocessing import cpu_count
import time

try:
    import requests
except:
    print("ModuleNotFoundError: No module named 'requests' (module 'requests' is not installed)")
    print("You can try install it by command:")
    print("pip3 install requests")
    print("or you can you can follow the steps described on:")
    print("https://onedata4sci.readthedocs.io/en/latest/user/onedata-downloader.html")
    sys.exit(1)

"""
Default verbosity level.
"""
VERBOSITY = 0

"""
Default Onezone service URL.
"""
DEFAULT_ONEZONE = "https://datahub.egi.eu"

"""
Used Onezone API URI.
"""
ONEZONE_API = "/api/v3/onezone/"

"""
Chunk size for downloading files as stream.
"""
CHUNK_SIZE = 33_554_432

"""
File extension of not yet completely downloaded (part) file.
"""
PART_FILE_EXTENSION = ".oddown_part"

"""
Threads number for parallel downloading.
"""
THREADS_NUMBER = 1

"""
Tries of downloading the file after error occurred
"""
TRIES_NUMBER = 3

"""
Number of seconds between two tries to download the file
"""
TRIES_DELAY = 5

ROOT_DIRECTORY_SIZE = 0

ALL_DIRECTORIES = 0
DIRECTORIES_CREATED = 0
DIRECTORIES_NOT_CREATED_OS_ERROR = 0

ALL_FILES = 0
EXISTENT_FILES = queue.Queue()
FINISHED_FILES = queue.Queue()
PART_FILES = queue.Queue()


file_queue = queue.Queue()
ERROR_QUEUE = queue.Queue()


def convert_chunk_size(chunk_size: str) -> int:
    """
    Converts user-given chunk size to integer.
    User can input values as number (bytes) or number + unit (eg. 32M)
    """
    chunk_size = chunk_size.strip()
    unit = "b"
    if chunk_size[-1].isalpha():
        unit = chunk_size[-1]
        chunk_size = chunk_size[:-1]

    try:
        chunk_size = int(chunk_size)
    except ValueError as e:
        print(
            "failed while converting size to integer, exception occured:",
            e.__class__.__name__,
        )
        return -1

    units = ("b", "k", "M", "G")
    if unit not in units:
        print("failed while converting mapping unit, unit is not in the right format")
        return -1

    unit_power = units.index(unit)
    chunk_size = chunk_size * (1024**unit_power)

    return chunk_size


def generate_random_string(size: int = 16) -> str:
    """
    Generates random string of characters of given size
    """
    if size < 0:
        return ""

    characters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890"
    random_string = "".join(random.choices(characters, k=size))
    return random_string


def remove_part_files(directory_to_search: str) -> bool:
    """
    Removes files in tree with extension defined by global value PART_FILE_EXTENSION
    """
    pattern = ".*\\" + PART_FILE_EXTENSION + "$"
    try:
        for root, directories, files in os.walk(directory_to_search):
            for actual_file in files:
                if re.match(pattern, actual_file):
                    os.remove(os.path.join(root, actual_file))
                    print(
                        "Partially downloaded file",
                        root + os.sep + actual_file,
                        "removed",
                    )
    except Exception as e:
        print("failed while removing part files, exception occured:", e.__class__.__name__)
        return False

    return True


def verbose_print(level, *args, **kwargs):
    """
    Print only when VERBOSITY is equal or higher than given level.
    """
    if VERBOSITY >= level:
        print(*args, **kwargs)


def error_printer(response: requests.Response, thread_number: int, directory: str, file_name: str):
    print("failed", end="")
    response_json = response.json()
    if (
            "error" in response_json
            and "details" in response_json["error"]
            and "errno" in response_json["error"]["details"]
            and "eacces" in response_json["error"]["details"]["errno"]
    ):
        verbose_print(1, f"Thread {thread_number}:", end=" ")
        print(
            "Downloading of",
            directory + os.sep + file_name,
            "failed, response error: permission denied",
        )
    elif (
            "error" in response_json
            and "details" in response_json["error"]
            and "errno" in response_json["error"]["details"]
            and "enoent" in response_json["error"]["details"]["errno"]
    ):
        verbose_print(1, f"Thread {thread_number}:", end=" ")
        print(
            "Downloading of",
            directory + os.sep + file_name,
            "failed, response error: no such file or directory",
        )
    else:
        verbose_print(1, f"Thread {thread_number}:", end=" ")
        print(
            "Downloading of",
            directory + os.sep + file_name,
            "failed, returned HTTP response code =",
            response.status_code,
        )

    verbose_print(1, f"Thread {thread_number}:", response.json())


def chunkwise_downloader(request: requests.Response, directory: str, file_name: str, random_filename: str, thread_number: int) -> int:
    try:
        with open(directory + os.sep + random_filename, "wb") as file:
            for chunk in request.iter_content(chunk_size=CHUNK_SIZE):
                file.write(chunk)
        os.rename(
            directory + os.sep + random_filename,
            directory + os.sep + file_name,
        )
        FINISHED_FILES.put(os.path.join(directory, file_name))
        verbose_print(
            2,
            f"Thread {thread_number}",
            random_filename,
            "renamed to",
            directory + os.sep + file_name,
        )
        verbose_print(1, f"Thread {thread_number}:", end=" ")
        print(
            f"Downloading file",
            directory + os.sep + file_name,
            "was successful",
        )
        return 0
    except EnvironmentError as e:
        verbose_print(1, f"Thread {thread_number}:", end=" ")
        print(f"Failed {directory + os.sep + file_name}, exception occured:", e.__class__.__name__)
        verbose_print(1, str(e))
        return 1


def download_file(onezone, file_id, file_name, directory, thread_number):
    """
    Download file with given file_id to given directory.
    """
    verbose_print(2, "download_file(%s, %s, %s, %s)" % (onezone, file_id, file_name, directory))
    # don't download the file when it exists
    random_filename = generate_random_string(size=16) + PART_FILE_EXTENSION
    verbose_print(1, f"Thread {thread_number}:", end=" ")
    verbose_print(1, "Downloading file", directory + os.sep + file_name, end=" ")
    verbose_print(2, " (temporary filename " + random_filename + ") ", end="")
    verbose_print(1, "started", flush=True)

    if os.path.exists(directory + os.sep + file_name):
        EXISTENT_FILES.put(os.path.join(directory, file_name))
        verbose_print(1, f"Thread {thread_number}:", end=" ")
        print("File", directory + os.sep + file_name, "exists, skipped")
        return 0

    # https://onedata.org/#/home/api/stable/oneprovider?anchor=operation/download_file_content
    url = onezone + ONEZONE_API + "shares/data/" + file_id + "/content"

    try_number = 0
    for try_number in range(TRIES_NUMBER):
        with requests.get(url, allow_redirects=True, stream=True) as request:
            if not request.ok:
                error_printer(request, thread_number, directory, file_name)
                return 2

            if try_number == 0:
                PART_FILES.put(os.path.join(directory, random_filename))

            if chunkwise_downloader(request, directory, file_name, random_filename, thread_number) == 0:
                break

            time.sleep(TRIES_DELAY)

    if try_number == TRIES_NUMBER - 1:
        ERROR_QUEUE.put(f"The file {directory + os.sep + file_name} could not be downloaded")
        return 2

    return 0


def process_directory(onezone, file_id, file_name, directory):
    """
    Process directory and recursively its content.
    """
    verbose_print(
        2,
        "process_directory(%s, %s, %s, %s)" % (onezone, file_id, file_name, directory),
    )
    global ALL_DIRECTORIES
    global DIRECTORIES_CREATED
    global DIRECTORIES_NOT_CREATED_OS_ERROR
    # don't create the the directory when it exists
    print("Processing directory", directory + os.sep + file_name, end="... ", flush=True)
    try:
        os.mkdir(directory + os.sep + file_name, mode=0o777)
        DIRECTORIES_CREATED += 1
        verbose_print(1, "directory created", end="")
        print()
    except FileExistsError:
        print("directory exists, not created", end="")
        print()
    except FileNotFoundError as e:
        DIRECTORIES_NOT_CREATED_OS_ERROR += 1
        print("failed, exception occured:", e.__class__.__name__)
        verbose_print(1, str(e))
        return 2

    # get content of new directory
    # https://onedata.org/#/home/api/stable/oneprovider?anchor=operation/list_children
    url = onezone + ONEZONE_API + "shares/data/" + file_id + "/children"
    response = requests.get(url)
    if response.ok:
        response_json = response.json()
        result = 0
        # process child nodes
        for child in response_json["children"]:
            # difference between Onezone version 20 and 21 in name of the key containing the file_id attribute
            if "file_id" in child:
                child_file_id = child["file_id"]
            else:
                child_file_id = child["id"]
            result = process_node(onezone, child_file_id, directory + os.sep + file_name) or result

        return result
    else:
        print("Error: failed to process directory", file_name)
        verbose_print(1, "processed directory", file_name, " with File ID =", file_id)
        verbose_print(1, response.json())
        return 2


def process_node(onezone, file_id, directory):
    """
    Process given node (directory or file).
    """
    verbose_print(2, "process_node(%s, %s, %s)" % (onezone, file_id, directory))
    global ROOT_DIRECTORY_SIZE
    global ALL_FILES
    global ALL_DIRECTORIES
    # get basic node's attributes
    # https://onedata.org/#/home/api/stable/oneprovider?anchor=operation/get_attrs
    url = onezone + ONEZONE_API + "shares/data/" + file_id
    response = requests.get(url)
    if response.ok:
        response_json = response.json()
        node_type = response_json["type"].upper()
        node_name = response_json["name"]
        node_size = response_json["size"]

        if node_size > ROOT_DIRECTORY_SIZE:
            ROOT_DIRECTORY_SIZE = node_size

        result = 0
        # check if node is directory or folder
        if node_type == "REG" or node_type == "SYMLNK":
            ALL_FILES += 1
            if os.path.exists(os.path.join(directory, node_name)):
                EXISTENT_FILES.put(os.path.join(directory, node_name))
                print(
                    "File",
                    directory + os.sep + node_name,
                    "exists, will not be downloaded",
                )
                return 0

            verbose_print(1, "Adding file to queue", directory + os.sep + node_name)
            file_queue.put((onezone, file_id, node_name, directory))
        elif node_type == "DIR":
            ALL_DIRECTORIES += 1
            result = process_directory(onezone, file_id, node_name, directory) or result
        else:
            print("Error: unknown node type")
            verbose_print(1, "returned node type", node_type, " of node with File ID =", file_id)
            verbose_print(1, response.json())
            return 2

        return result
    else:
        print(
            "Error: failed to retrieve information about the node. The requested node may not exist."
        )
        verbose_print(1, "requested node File ID =", file_id)
        verbose_print(1, response.json())
        return 1


def clean_onezone(onezone):
    """
    Clean and test of given Onezone service.
    """
    # add protocol name if not specified
    if not onezone.startswith("https://") and not onezone.startswith("http://"):
        onezone = "http://" + onezone

    verbose_print(1, "Use Onezone:", onezone)

    # test if such Onezone exists
    url = onezone + ONEZONE_API + "configuration"
    try:
        response = requests.get(url)
    except Exception as e:
        print("Error: failure while trying to communicate with Onezone:", onezone)
        verbose_print(1, str(e))
        sys.exit(2)

    if not response.ok:
        print("Error: failure while connecting to Onezone:", onezone)
        sys.exit(2)

    try:
        response_json = response.json()
        verbose_print(2, "Onezone configuration:")
        verbose_print(2, response_json)
    except Exception as e:
        print(
            "Error: failure while parsing JSON response from Onezone:",
            e.__class__.__name__,
        )
        sys.exit(2)

    # get Onezone version
    onezone_version = response_json["version"].split(".")[0]  # 21.02.0-alpha28 -> 21
    verbose_print(1, "Onezone version:", onezone_version)

    return onezone


def clean_directory(directory):
    """
    Test if given directory is correct.
    """
    # test if given directory exists
    if not os.path.isdir(directory):
        print("Error: output directory", directory, "does not exist")
        sys.exit(2)

    return directory


def thread_worker(thread_number: int):
    result = 0
    while True:
        onezone, file_id, node_name, directory = file_queue.get()
        result = download_file(onezone, file_id, node_name, directory, thread_number) or result
        file_queue.task_done()
    return result


def print_download_statistics(directory_to_search: str, finished: bool = True):
    errors = ERROR_QUEUE.qsize()

    existent_files = EXISTENT_FILES.qsize()
    finished_files = FINISHED_FILES.qsize()

    part_size = 0
    while not PART_FILES.empty():
        file_path = PART_FILES.get()
        if os.path.exists(file_path):
            part_size += os.path.getsize(file_path)

    finished_size = 0
    while not FINISHED_FILES.empty():
        file_path = FINISHED_FILES.get()
        finished_size += os.path.getsize(file_path)

    existent_size = 0
    while not EXISTENT_FILES.empty():
        file_path = EXISTENT_FILES.get()
        existent_size += os.path.getsize(file_path)

    downloaded_size = finished_size + part_size

    print()
    if errors != 0:
        print("Errors during execution:")
        while not ERROR_QUEUE.empty():
            print(ERROR_QUEUE.get())
        print()

    print("Download statistics:")
    if ALL_FILES != 0:
        print(
            f"Files created: {finished_files}/{ALL_FILES} ({(finished_files/ALL_FILES * 100):.2f}%), already existent: {existent_files}, error while creating: {ALL_FILES - (existent_files + finished_files)}"
        )
    else:
        print(
            f"Files created: 0, already existent: {existent_files}, error while creating: {ALL_FILES - (existent_files + finished_files)}"
        )
    if ALL_DIRECTORIES != 0:
        print(
            f"Directories created: {DIRECTORIES_CREATED}/{ALL_DIRECTORIES} ({DIRECTORIES_CREATED/ALL_DIRECTORIES * 100:.2f}%), already existent: {ALL_DIRECTORIES - (DIRECTORIES_NOT_CREATED_OS_ERROR + DIRECTORIES_CREATED)}, error while creating: {DIRECTORIES_NOT_CREATED_OS_ERROR}"
        )
    else:
        print(
            f"Directories created: 0, already existent: {ALL_DIRECTORIES - (DIRECTORIES_NOT_CREATED_OS_ERROR + DIRECTORIES_CREATED)}, error while creating: {DIRECTORIES_NOT_CREATED_OS_ERROR}"
        )
    if ROOT_DIRECTORY_SIZE != 0:
        print(
            f"Downloaded size: {downloaded_size}/{ROOT_DIRECTORY_SIZE} bytes ({downloaded_size/ROOT_DIRECTORY_SIZE * 100:.2f}%), finished: {finished_size} bytes, existent: {existent_size} bytes, part files: {part_size} bytes, not downloaded yet or error: {ROOT_DIRECTORY_SIZE - (finished_size + existent_size + part_size)} bytes"
        )
    else:
        print(
            f"Downloaded size: 0 bytes, finished: {finished_size} bytes, existent: {existent_size} bytes, part files: {part_size} bytes, not downloaded yet or error: {ROOT_DIRECTORY_SIZE - (finished_size + existent_size + part_size)} bytes"
        )
    if not finished:
        print("RESULTS MAY BE INCORRECT, PROGRAM DID NOT FINISH CORRECTLY")


def main():
    parser = argparse.ArgumentParser(
        description="Script allowing you to download a complete shared space, directory or even a single file from the Onedata system."
    )
    parser.add_argument(
        "-o",
        "--onezone",
        default=DEFAULT_ONEZONE,
        type=str,
        help="Onedata Onezone URL with specified protocol (default: https://datahub.egi.eu)",
    )
    parser.add_argument(
        "-d",
        "--directory",
        default=".",
        type=str,
        help="Output directory (default: current directory)",
    )
    parser.add_argument(
        "-c",
        "--chunk-size",
        default="32M",
        type=str,
        help="The size of downloaded file segments (chunks) after which the file is written to disk. Value can be in bytes, or a number with unit e.g. 16k, 32M or 2G (default: 32M).",
    )
    parser.add_argument(
        "-j",
        "--threads-number",
        default=cpu_count(),
        type=int,
        help="Number of threads for parallel downloading. Setting this parameter to a reasonable value can significantly reduce the overall download time (default: 1).",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Set verbose prints - displaying debug information",
    )
    parser.add_argument("file_id", type=str, help="File ID of shared space, directory or a file")
    args = parser.parse_args()

    # set up verbosity level
    global VERBOSITY
    VERBOSITY = args.verbose

    global CHUNK_SIZE
    CHUNK_SIZE = convert_chunk_size(args.chunk_size)
    if CHUNK_SIZE < 0:
        return 3

    status = remove_part_files(args.directory)
    if not status:
        return 4

    global THREADS_NUMBER
    THREADS_NUMBER = args.threads_number
    if THREADS_NUMBER < 1:
        print("failed on startup; number of threads cannot be lower than one")
        return 5

    onezone = clean_onezone(args.onezone)
    directory = clean_directory(args.directory)

    try:
        print("Exploring and creating the directory structure")
        result = process_node(onezone, args.file_id, directory)
        if result:
            print_download_statistics(args.directory, finished=False)
            return result

        print("Downloading files")
        for thread_number in range(THREADS_NUMBER):
            result = (
                threading.Thread(target=thread_worker, args=(thread_number,), daemon=True).start()
                or result
            )
        file_queue.join()
        result = 0 if ERROR_QUEUE.qsize() == 0 else 1
        print_download_statistics(args.directory)
        return result
    except KeyboardInterrupt as e:
        print(" prematurely interrupted (" + e.__class__.__name__ + ")")
        print_download_statistics(args.directory, finished=False)
        return 2


if __name__ == "__main__":
    return_code = main()
    sys.exit(return_code)