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
from typing import Optional, Generator

try:
    import requests
except:
    print("ModuleNotFoundError: No module named 'requests' (module 'requests' is not installed)")
    print("You can try install it by command:")
    print("pip3 install requests")
    print("or you can you can follow the steps described on:")
    print("https://onedata4sci.readthedocs.io/en/latest/user/onedata-downloader.html")
    sys.exit(1)


class VERBOSITY:
    DEF = 0
    V = 1
    VV = 2


V = VERBOSITY

"""Max (lowest) priority for downloading the file 
"""
MAX_PRIORITY = 3

"""Tries of downloading the file after error occurred
"""
TRIES_NUMBER: int = 7

"""
Default verbosity level.
"""
VERBOSITY = V.DEF

"""
Default Onezone service URL.
"""
DEFAULT_ONEZONE = "https://datahub.egi.eu"

"""
Used Onezone API URI.
"""
ONEZONE_API: str = "/api/v3/onezone/"

"""
Chunk size for downloading files as stream in bytes.
"""
CHUNK_SIZE: int = 32 * 1024 * 1024  # 32 MB - 33_554_432

"""
File extension of not yet completely downloaded (part) file.
"""
PART_FILE_EXTENSION: str = ".oddown_part"

"""
Threads number for parallel downloading.
"""
THREADS_NUMBER: int = 1

"""
Number of seconds between two tries to download the file
"""
TRIES_DELAY: int = 1

ONEZONE: str = DEFAULT_ONEZONE

DIRECTORY: str = "."

FILE_ID: Optional[str] = None


def priority_subtractor():
    iterator = 0
    while True:
        yield -1 if (iterator == 0) else 0
        iterator = (iterator + 1) % 3


class URLs:
    def __init__(self, onezone: str, file_id: str):
        self._content = onezone + ONEZONE_API + "shares/data/" + file_id + "/content"
        self._children = onezone + ONEZONE_API + "shares/data/" + file_id + "/children"
        self._node_attributes = onezone + ONEZONE_API + "shares/data/" + file_id

    @property
    def content(self):
        # https://onedata.org/#/home/api/stable/oneprovider?anchor=operation/download_file_content
        return self._content

    @property
    def children(self):
        # https://onedata.org/#/home/api/stable/oneprovider?anchor=operation/list_children
        return self._children

    @property
    def node_attrs(self):
        # https://onedata.org/#/home/api/stable/oneprovider?anchor=operation/get_attrs
        return self._node_attributes


class DownloadableItem(object):
    def __init__(self, onezone: str, file_id: str, node_name: str, directory: str):
        self._onezone: str = onezone
        self._file_id: str = file_id
        self._node_name: str = node_name
        self._directory: str = directory
        self._priority: int = MAX_PRIORITY  # internal value, lowering
        self._ttl: int = TRIES_NUMBER
        self._part_filename: str = generate_random_string(size=16) + PART_FILE_EXTENSION
        self._priority_subtractor: Generator = priority_subtractor()
        self._path = os.path.join(self._directory, self._node_name)  # not to compute it again
        self._part_path = os.path.join(self._directory, self._part_filename)  # not to compute it again
        self._urls = URLs(self._onezone, self._file_id)

    @property
    def onezone(self) -> str:
        return self._onezone

    @property
    def file_id(self) -> str:
        return self._file_id

    @property
    def node_name(self) -> str:
        return self._node_name

    @property
    def directory(self) -> str:
        return self._directory

    @property
    def path(self) -> str:
        return self._path

    @property
    def priority(self) -> int:
        """Number representing priority, lower number is higher priority
        """
        return MAX_PRIORITY - self._priority

    @property
    def part_filename(self) -> str:
        return self._part_filename

    @property
    def part_path(self) -> str:
        return self._part_path

    @property
    def URL(self) -> URLs:
        return self._urls

    def _decrease_priority(self) -> None:
        """Lowers the priority by one step
        """
        self._priority = max(0, self._priority + next(self._priority_subtractor))

    def try_to_download(self) -> bool:
        if self._ttl == 0:
            return False

        self._ttl -= 1
        self._decrease_priority()
        return True

    def __lt__(self, other) -> bool:
        """Sorting by priority, lower number is higher priority
        """
        if not isinstance(other, DownloadableItem):
            raise TypeError(f"Instance of {type(other).__name__} cannot be compared with {type(self).__name__}")
        return self.priority < other.priority


class QueuePool:
    def __init__(self, queues: tuple[queue.Queue, ...], weights: tuple[int, ...]):
        if len(queues) != len(weights):
            raise AttributeError("Number of queues must be equal to number of their weights")

        self._queues = queues
        self._weights = weights

        _weights = []
        for item in [[index] * weight for index, weight in enumerate(weights)]:
            _weights.extend(item)
        random.shuffle(_weights)  # dequeue to be more fair

        self._weight_queue = queue.Queue()
        for weight in _weights:
            self._weight_queue.put(weight)

        self._queue_to_finish = 0
        self._mutex = threading.Lock()

    def __len__(self):
        return len(self._queues)

    def join(self):
        for key, act_queue in enumerate(self._queues):
            act_queue.join()

    def _increase_weight(self, index: int, thread_number: int):
        if index == len(self) - 1:
            return

        passed_queue_items = 0
        while (self._weight_queue.qsize() != sum(self._weights[self._queue_to_finish + 1:])
               and passed_queue_items < sum(self._weights)):
            # deadlock possibility if different thread has not yet put number back
            try:
                act_ind = self._weight_queue.get(block=False)
            except queue.Empty:
                passed_queue_items += 1
                continue
            if act_ind != index:
                self._weight_queue.put(act_ind)
            passed_queue_items += 1

        if self._queue_to_finish == index:
            v_print(V.V, f"Thread {thread_number}: increasing queue index {self._queue_to_finish}++")
            self._queue_to_finish += 1

    def _try_to_increase_weight(self, index: int, thread_number: int) -> int:
        v_print(V.VV, f"Thread {thread_number}: trying to acquire mutex")
        self._mutex.acquire(blocking=True)
        v_print(V.VV, f"Thread {thread_number}: mutex acquired")

        if not (self._queue_to_finish >= index and self.get_queue(index).qsize() == 0):
            v_print(V.VV, f"Thread {thread_number}: condition not met, releasing")
            self._mutex.release()
            index = self._weight_queue.get()
            return index

        self._increase_weight(index, thread_number)

        v_print(V.VV, f"Thread {thread_number}: releasing mutex")
        self._mutex.release()

        index = self._weight_queue.get()
        return index

    def fair_index(self, thread_number: int):
        index = self._weight_queue.get()

        if self._queue_to_finish >= index and self.get_queue(index).qsize() == 0:
            self._weight_queue.put(index)
            index = self._try_to_increase_weight(index, thread_number)

        self._weight_queue.put(index)
        return index

    def get_queue(self, index: int) -> queue.Queue:
        if index >= len(self._queues) or index < 0:
            raise IndexError("Queue index out of bounds")
        return self._queues[index]


ROOT_DIRECTORY_SIZE = 0
ALL_DIRECTORIES = 0
DIRECTORIES_CREATED = 0
DIRECTORIES_NOT_CREATED_OS_ERROR = 0

ALL_FILES = 0
EXISTENT_FILES = queue.Queue()
FINISHED_FILES = queue.Queue()
PART_FILES = queue.Queue()


_file_queue = queue.Queue()
_priority_file_queue = queue.PriorityQueue()
QP = QueuePool(queues=(_file_queue, _priority_file_queue), weights=(15, 1))

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
        v_print(V.DEF, f"failed while converting size to integer, exception occured: {e.__class__.__name__}")
        return -1

    units = ("b", "k", "M", "G")
    if unit not in units:
        v_print(V.DEF, "failed while converting mapping unit, unit is not in the right format")
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
                    file_path = os.path.join(root, actual_file)
                    try:
                        os.remove(file_path)  # cannot get OSError, because not going through directories
                    except FileNotFoundError:
                        v_print(V.DEF, f"cannot remove {file_path}, it does not exist")
                    else:
                        v_print(V.DEF, f"Partially downloaded file {file_path} removed")
    except OSError as e:
        v_print(V.DEF, "failed while removing part files, exception occured:", e.__class__.__name__)
        return False

    return True


def verbose_print(level, *args, **kwargs):
    """
    Print only when VERBOSITY is equal or higher than given level.
    """
    if VERBOSITY >= level:
        print(*args, **kwargs)


v_print = verbose_print


def error_printer(response: requests.Response, thread_number: int, file: DownloadableItem):
    v_print(V.DEF, "failed", end="")
    response_json = response.json()

    v_print(V.V, f"Thread {thread_number}:", end=" ")
    if (
            "error" in response_json  # dict
            and "details" in response_json["error"]  # dict
            and "errno" in response_json["error"]["details"]
            and response_json["error"]["details"]["errno"] in ("eaccess", "enoent")
    ):
        if "eacces" in response_json["error"]["details"]["errno"]:
            v_print(V.DEF, f"Downloading of {file.path} failed, response error: permission denied")
        if "enoent" in response_json["error"]["details"]["errno"]:
            v_print(V.DEF, f"Downloading of {file.path} failed, response error: no such file or directory")
    else:
        v_print(V.DEF, f"Downloading of {file.path} failed, returned HTTP response code = {response.status_code}")

    v_print(V.V, f"Thread {thread_number}:", response_json)


def chunkwise_downloader(request: requests.Response, file: DownloadableItem, thread_number: int) -> int:
    try:
        with open(file.part_path, "ab") as f:  # if file was already opened and written into, it will continue
            for chunk in request.iter_content(chunk_size=CHUNK_SIZE, decode_unicode=True):
            # for chunk in request.iter_content():
                f.write(chunk)
                # flushing automatically as OS says
        # the file is closed now
    except EnvironmentError as e:
        v_print(V.V, f"Thread {thread_number}:", end=" ")
        v_print(V.DEF, f"Failed {file.path}, exception occured:", e.__class__.__name__)
        v_print(V.V, str(e))
        return 1

    return 0


def renamer(file: DownloadableItem, thread_number: int):
    try:
        os.rename(file.part_path, file.path)
        FINISHED_FILES.put(file.path)

        v_print(V.VV, f"Thread {thread_number}: {file.part_filename} renamed to {file.path}")
        v_print(V.V, f"Thread {thread_number}:", end=" ")
    except OSError:
        v_print(V.V, f"Thread {thread_number}: could not rename {file.part_filename} to {file.path}")
        return 1

    return 0


def download_file(file: DownloadableItem, thread_number: int):
    """
    Download file with given file_id to given directory.
    """
    v_print(V.VV, f"download_file({file.onezone}, {file.file_id}, {file.node_name}, {file.directory})")
    # don't download the file when it exists

    v_print(V.V, f"Thread {thread_number}:", end=" ")
    v_print(V.V, "Downloading file", file.path, end=" ")
    v_print(V.VV, " (temporary filename " + file.part_filename + ") ", end="")
    v_print(V.V, "started", flush=True)

    if os.path.exists(file.path):
        EXISTENT_FILES.put(file.path)
        v_print(V.V, f"Thread {thread_number}:", end=" ")
        v_print(V.DEF, "File", file.path, "exists, skipped")
        return 0

    headers = {}
    already_downloaded = 0
    if os.path.exists(file.part_path):  # incorrectly downloaded
        v_print(V.VV, f"Thread {thread_number}:", end=" ")
        v_print(V.V, f"part file exists ({file.part_path})", end=", ")
        already_downloaded = os.path.getsize(file.part_path)
        v_print(V.V, f"already downloaded {already_downloaded} bytes")
        headers["Range"] = f"bytes={already_downloaded}-"

    with requests.get(file.URL.content, headers=headers, allow_redirects=True, stream=True) as request:
        if request.status_code == 416:
            v_print(V.VV, f"Thread {thread_number}:", end=" ")
            v_print(V.V, "got status code 416 while downloading, trying to get the original size")
            with requests.get(file.URL.content, allow_redirects=True, stream=True) as request_size:
                original_size = request_size.headers.get("content-length")
                if already_downloaded != original_size:
                    v_print(V.V, f"the original size does not match, already downloaded: {already_downloaded}, "
                                 f"file size: {original_size}")
                    return 5
                v_print(V.V, f"the original size does matches, the size is: {already_downloaded}")
        else:
            if not request.ok:
                error_printer(request, thread_number, file)
                return 2

            if chunkwise_downloader(request, file, thread_number) != 0:
                return 3

    if renamer(file, thread_number) != 0:
        return 4

    v_print(V.DEF, f"Downloading file {file.path} was successful")

    return 0


def process_directory(onezone, file_id, file_name, directory):
    """
    Process directory and recursively its content.
    """
    v_print(V.VV,f"process_directory({onezone}, {file_id}, {file_name}, {directory})")
    global ALL_DIRECTORIES
    global DIRECTORIES_CREATED
    global DIRECTORIES_NOT_CREATED_OS_ERROR
    # don't create the the directory when it exists
    v_print(V.DEF, "Processing directory", directory + os.sep + file_name, end="... ", flush=True)
    try:
        os.mkdir(directory + os.sep + file_name, mode=0o777)
        DIRECTORIES_CREATED += 1
        v_print(V.V, "directory created")
    except FileExistsError:  # directory already existent
        v_print(V.DEF, "directory exists, not created")
    except FileNotFoundError as e:  # parent directory non existent
        DIRECTORIES_NOT_CREATED_OS_ERROR += 1
        v_print(V.DEF, "failed, exception occured:", e.__class__.__name__)
        v_print(V.V, str(e))
        return 2

    # get content of new directory

    response = requests.get(URLs(onezone, file_id).children)
    if not response.ok:
        v_print(V.DEF, "Error: failed to process directory", file_name)
        v_print(V.V, "processed directory", file_name, " with File ID =", file_id)
        v_print(V.V, response.json())
        return 2

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


def process_node(onezone: str, file_id: str, directory: str):
    """
    Process given node (directory or file).
    """
    v_print(V.VV, "process_node(%s, %s, %s)" % (onezone, file_id, directory))
    global ROOT_DIRECTORY_SIZE
    global ALL_FILES
    global ALL_DIRECTORIES
    # get basic node's attributes

    response = requests.get(URLs(onezone, file_id).node_attrs)
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
            node_path = os.path.join(directory, node_name)
            if os.path.exists(node_path):
                EXISTENT_FILES.put(node_path)
                v_print(V.DEF, f"File {node_path} exists, it will not be downloaded")
                return 0

            v_print(V.V, "Adding file to queue", node_path)
            file_queue = QP.get_queue(0)
            file_queue.put(DownloadableItem(onezone, file_id, node_name, directory))
        elif node_type == "DIR":
            ALL_DIRECTORIES += 1
            result = process_directory(onezone, file_id, node_name, directory) or result
        else:
            v_print(V.DEF, "Error: unknown node type")
            v_print(V.V, "returned node type", node_type, " of node with File ID =", file_id)
            v_print(V.V, response.json())
            return 2

        return result
    else:
        v_print(V.DEF,
            "Error: failed to retrieve information about the node. The requested node may not exist."
        )
        v_print(V.V, "requested node File ID =", file_id)
        v_print(V.V, response.json())
        return 1


def clean_onezone(onezone):
    """
    Clean and test of given Onezone service.
    """
    # add protocol name if not specified
    if not onezone.startswith("https://") and not onezone.startswith("http://"):
        onezone = "http://" + onezone

    v_print(V.V, "Use Onezone:", onezone)

    # test if such Onezone exists
    url = onezone + ONEZONE_API + "configuration"
    try:
        response = requests.get(url)
    except Exception as e:
        v_print(V.DEF, "Error: failure while trying to communicate with Onezone:", onezone)
        v_print(V.V, str(e))
        sys.exit(2)

    if not response.ok:
        v_print(V.DEF, "Error: failure while connecting to Onezone:", onezone)
        sys.exit(2)

    try:
        response_json = response.json()
        v_print(V.VV, "Onezone configuration:")
        v_print(V.VV, response_json)
    except Exception as e:
        v_print(V.DEF,
            "Error: failure while parsing JSON response from Onezone:",
            e.__class__.__name__,
        )
        sys.exit(2)

    # get Onezone version
    onezone_version = response_json["version"].split(".")[0]  # 21.02.0-alpha28 -> 21
    v_print(V.V, "Onezone version:", onezone_version)

    return onezone


def clean_directory(directory):
    """
    Test if given directory is correct.
    """
    # test if given directory exists
    if not os.path.isdir(directory):
        v_print(V.DEF, "Error: output directory", directory, "does not exist")
        sys.exit(2)

    return directory


def thread_worker(thread_number: int):
    while True:
        queue_index = QP.fair_index(thread_number)
        actual_queue = QP.get_queue(queue_index)
        try:
            v_print(V.V, f"Thread {thread_number}: acquiring download or blocked state in queue {queue_index}")
            downloadable_item: DownloadableItem = actual_queue.get(block=True)
            v_print(V.V, f"Thread {thread_number}: acquired download in {queue_index}")
        except queue.Empty:  # incorrect queue
            continue

        if queue_index == 0:
            PART_FILES.put(downloadable_item.part_path)

        v_print(V.VV, f"Thread: {thread_number}, actual queue index: {queue_index}, file priority: {downloadable_item.priority}, ttl: {downloadable_item._ttl}")
        if downloadable_item.try_to_download():
            result = download_file(downloadable_item, thread_number)

            if result != 0:
                QP.get_queue(1).put(downloadable_item)
        else:
            ERROR_QUEUE.put(f"The file {downloadable_item.path} could not be downloaded")

        actual_queue.task_done()

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


def setup_parser() -> argparse.ArgumentParser:
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
        default=1,
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

    return parser


def process_parser(parser: argparse.ArgumentParser):
    args = parser.parse_args()

    # set up verbosity level
    global VERBOSITY
    VERBOSITY = args.verbose

    global CHUNK_SIZE
    CHUNK_SIZE = convert_chunk_size(args.chunk_size)
    if CHUNK_SIZE < 0:
        return 3

    global THREADS_NUMBER
    THREADS_NUMBER = args.threads_number
    if THREADS_NUMBER < 1:
        v_print(V.DEF, "failed on startup; number of threads cannot be lower than one")
        return 4

    global ONEZONE
    ONEZONE = clean_onezone(args.onezone)

    global DIRECTORY
    DIRECTORY = clean_directory(args.directory)

    global FILE_ID
    FILE_ID = args.file_id

    return 0


def main():
    parser = setup_parser()
    result = process_parser(parser)
    if result != 0:
        return result

    result = remove_part_files(DIRECTORY)
    if not result:
        return 5

    try:
        v_print(V.DEF, "Exploring and creating the directory structure")
        result = process_node(ONEZONE, FILE_ID, DIRECTORY)
        if result:
            print_download_statistics(DIRECTORY, finished=False)
            return result

        v_print(V.DEF, "Downloading files")
        for thread_number in range(THREADS_NUMBER):
            result = (
                threading.Thread(target=thread_worker, args=(thread_number,), daemon=True).start()
                or result
            )
        QP.join()
        result = 0 if ERROR_QUEUE.qsize() == 0 else 1
        print_download_statistics(DIRECTORY)
        return result
    except KeyboardInterrupt as e:
        v_print(V.DEF, " prematurely interrupted (" + e.__class__.__name__ + ")")
        print_download_statistics(DIRECTORY, finished=False)
        return 2


if __name__ == "__main__":
    return_code = main()
    sys.exit(return_code)
