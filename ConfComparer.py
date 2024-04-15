#!/usr/bin/python
import argparse
import csv
import sys
from enum import Enum, auto
import mmap
from pathlib import Path
from shutil import which
import subprocess
import tempfile
from typing import Callable, Dict, List, Optional
import warnings

##############################
# CONSTANTS AND CUSTOM TYPES #
##############################
# supported extensions in this list must be UPPERCASE
_SUPPORTED_EXTENSIONS: List[str] = [".PDB"]


##############
# EXCEPTIONS #
##############
class ConfComparerError(Exception):
    """General exception for this project"""


################
# HELPER STUFF #
################
class MoleculeNameCase(Enum):
    NO_CHANGE = auto()
    LOWER = auto()
    UPPER = auto()


_CASE_TRANSFORMING_FUNCTIONS: Dict[MoleculeNameCase, Callable[[str], str]] = {
    MoleculeNameCase.NO_CHANGE: (lambda s: s),
    MoleculeNameCase.LOWER: str.lower,
    MoleculeNameCase.UPPER: str.upper
}


########################
# RUNTIME DATA STORAGE #
########################
class _Storage:
    _script_folder = Path(__file__).resolve().parent
    templates_dir: Path = _script_folder.joinpath("templates")
    input_dir: Path = _script_folder.joinpath("input")
    output_dir: Path = _script_folder.joinpath("output")
    SB_executable_path: Path = _script_folder.joinpath("SB_batch", "SiteBinderCMD.exe")
    use_mono: bool = False
    tolerance: float = 1.0
    print_list = False
    print_RMSD_chart = False
    print_summary = True


#################################
# ESSENTIAL CLASSES AND OBJECTS #
#################################
class MoleculeBase:
    """
    Represents any molecule, carries info about name and path to source file
    """
    def __init__(self,
                 molecule_file_path: Optional[Path],
                 name_case: MoleculeNameCase = MoleculeNameCase.NO_CHANGE):
        self._case_transforming_func = _CASE_TRANSFORMING_FUNCTIONS[name_case]
        self.path = self.extension = None
        self.name = "none"
        self.ligand = ""

        if molecule_file_path is None:
            return

        path = molecule_file_path.resolve()
        ext = path.suffix
        if path and ext.upper() not in _SUPPORTED_EXTENSIONS:
            raise ConfComparerError(f"file '{path}': '{ext}' is not a supported file extension")

        # SB uses file name without extension(s) as molecule name
        self.name = path.name[:-len(ext)]
        self.path = str(path)
        self.extension = ext

    def get_name(self) -> str:
        return self._case_transforming_func(self.name)

    @classmethod
    def get_ligand_from_pdb(cls, molecule_path: Path) -> str:
        with molecule_path.open("rb") as f:
            with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as s:
                atom_record_start = s.find(b"HETATM")
                if atom_record_start == -1:
                    atom_record_start = s.find(b"ATOM")
                if atom_record_start == -1:
                    return ""
                s.seek(atom_record_start + 17)
                return s.read(3).decode()


class Conformation(MoleculeBase):
    """Represents particular conformation, carries info about name and path to template."""
    _NoneConformation: Optional["Conformation"] = None

    def __init__(self, template_file_path: Optional[Path]):
        super().__init__(template_file_path, MoleculeNameCase.UPPER)

    @classmethod
    def get_none_conformation(cls) -> "Conformation":
        if cls._NoneConformation is None:
            cls._NoneConformation = cls(None)
        return cls._NoneConformation


class Molecule(MoleculeBase):
    """
    Represents particular molecule, carries info about name, path to source file and computed
    RMSD values for various templates.
    """
    def __init__(self, molecule_file_path: Optional[Path]):
        super().__init__(molecule_file_path)
        if self.extension.upper() == ".PDB" and molecule_file_path:
            self.ligand = self.get_ligand_from_pdb(molecule_file_path)
        self.conformation_map: Dict[str, float] = {}

    @property
    def conformation(self) -> str:
        try:
            conf, *_ = {conf: rmsd for conf, rmsd in sorted(self.conformation_map.items(),
                                                            key=lambda item: item[1])
                        if rmsd <= _Storage.tolerance}
            return conf
        except ValueError:
            return Conformation.get_none_conformation().get_name()

    def update_conformation_map(self, conformation_name: str, rmsd: float):
        self.conformation_map[conformation_name] = rmsd


##############
# CORE LOGIC #
##############
class Analysis:
    def __init__(self):
        self._prepare_environment()

        self.conformations: Dict[str, Conformation] = dict()
        self.molecules: Dict[str, Molecule] = dict()

        self._expected_result_header = '"File","Rmsd","MatchedCount","ASize","BSize"'

        self.result_list_path = _Storage.output_dir.joinpath("results_individual_molecules.txt")
        self.result_rmsd_chart_path = _Storage.output_dir.joinpath("result_rmsd_chart.csv")
        self.result_summary_path = _Storage.output_dir.joinpath("result_summary.txt")

    def execute(self):
        # initialize conformations
        print("Loading templates...", file=sys.stderr)
        for f in _Storage.templates_dir.rglob("*"):
            if f.suffix.upper() not in _SUPPORTED_EXTENSIONS:
                continue
            conf = Conformation(f.resolve())
            self.conformations[conf.get_name()] = conf
        print("Templates loaded.", file=sys.stderr)

        # initialize molecules
        print("Loading molecules...", file=sys.stderr)
        for f in _Storage.input_dir.rglob("*"):
            if f.suffix.upper() not in _SUPPORTED_EXTENSIONS:
                continue
            mol = Molecule(f.resolve())
            self.molecules[mol.get_name()] = mol
        print("Molecules loaded.", file=sys.stderr)

        # work in a temporary folder where all the intermediate files can be stored during
        # the execution
        with tempfile.TemporaryDirectory() as tmp_dir:
            for conf_name, conf in self.conformations.items():
                print(f"Starting SB analysis of conformation {conf_name}...", file=sys.stderr)
                conf_folder = Path(tmp_dir).joinpath(conf_name)
                conf_folder.mkdir(mode=0o777)

                # create list of paths of molecules to be tested, add path of template on the first
                # line (it will be used as reference and discarded by SiteBinder)
                conf_input_file = conf_folder.joinpath(f"sb_input.list")
                print(f"Generating input file {conf_input_file}...", file=sys.stderr)
                with conf_input_file.open("w") as f:
                    for file_path in [conf.path] + [m.path for m in self.molecules.values()]:
                        f.write(f"{file_path}\n")

                # Usage: SiteBinderCMD input.txt rmsd.csv pairing.csv
                rmsd_output_file = conf_folder.joinpath(f"sb_output_rmsd.csv")
                pairing_output_file = conf_folder.joinpath(f"sb_output_pairing.csv")
                command_args = [arg for arg
                                in (_Storage.crossplatform_runner, _Storage.SB_executable_path,
                                    conf_input_file, rmsd_output_file, pairing_output_file)
                                if arg]

                print(f"Starting SiteBinder...", file=sys.stderr)
                sb_process = subprocess.Popen(command_args, stdout=subprocess.PIPE,
                                              stderr=subprocess.PIPE, shell=False)
                sb_process.wait()
                stdout, stderr = [stream.decode().strip() for stream in sb_process.communicate()]
                sb_process.stdout.close()

                if stdout or stderr or sb_process.returncode != 0:
                    print(f"error within processing conformation '{conf_name}', SiteBinder message:"
                          f"\nstdout: {stdout}\nstderr: {stderr}\n---", file=sys.stderr)
                    sys.exit(1)

                # process results for this conformation
                with rmsd_output_file.open("r") as f:
                    # first line is a header
                    if f.readline().strip() != self._expected_result_header:
                        raise ConfComparerError(f"Result file '{rmsd_output_file}' is not in "
                                                f"expected format")
                    # rest of the file are actual data, one entry per line
                    for entry in list(csv.reader(f)):
                        self._process_entry(conf_name, entry)
                print(f"Analysis of conformation {conf_name} finished.", file=sys.stderr)
            self._generate_result_files()
            print(file=sys.stderr)
        self._print_results()

    def _generate_result_files(self):
        # simple list of molecules and theirs conformations
        with self.result_list_path.open("w") as f:
            f.writelines(
                [f"[{molecule.ligand or 'UNKNOWN'}] {name}: {molecule.conformation}\n"
                 for name, molecule in self.molecules.items()]
            )

        # complex CSV table of RMSD values of every molecule with every conformation template
        with self.result_rmsd_chart_path.open("w") as f:
            # header with conformation names
            f.write(";;")
            for conf_name, _ in self.conformations.items():
                f.write(f"{conf_name};")
            f.write("\n")

            # individual molecules
            for name, molecule in self.molecules.items():
                f.write(f"{molecule.ligand};{name};")
                # use conf name from self.conformations instead of molecule.conformation map
                # because of header order
                for conf_name, _ in self.conformations.items():
                    f.write(f"{molecule.conformation_map.get(conf_name, 'NaN')};")
                f.write("\n")

        # summary
        # get number of occurrences for every conformation
        conf_absolute_numbers = {}
        for molecule in self.molecules.values():
            conf = molecule.conformation
            conf_absolute_numbers[conf] = conf_absolute_numbers.get(conf, 0) + 1

        # get total number of molecules
        total_num = len(self.molecules)

        # generate summary
        with self.result_summary_path.open("w") as f:
            for conf_name in list(self.conformations) + [
                    Conformation.get_none_conformation().get_name()]:
                f.write(f"{conf_name:<13} {conf_absolute_numbers.get(conf_name, 0): 4} "
                        f"{conf_absolute_numbers.get(conf_name, 0) * 100 /total_num:.2f}%\n")
            f.write(f"{'TOTAL':<13} {total_num: 4} 100%")

    def _print_results(self):
        if _Storage.print_list:
            print(self.result_list_path.read_text())
        if _Storage.print_RMSD_chart:
            print(self.result_rmsd_chart_path.read_text())
        if _Storage.print_summary:
            print(self.result_summary_path.read_text())

    def _process_entry(self, conf_name: str, entry: List[str]) -> None:
        if len(entry) == 5:
            name, rmsd, matched_count, a_size, b_size = entry
        elif len(entry) == 6:
            # SB uses simple comma as both CSV separator and decimal point, what is
            # incredibly confusing...
            name, rmsd_1, rmsd_2, matched_count, a_size, b_size = entry
            rmsd = f"{rmsd_1}.{rmsd_2}"
        else:
            warnings.warn(f"Skipping entry '{entry}' as the list contains unexpected number "
                          f"of members that cannot be matched to expected header "
                          f"'{self._expected_result_header}'")
            return
        if int(matched_count) != int(a_size) != int(b_size):
            warnings.warn(f"Molecule '{name}' is not of the same type as template "
                          f"of conformation '{conf_name}' (they probably differ in "
                          f"size or contain incompatible atoms.")
        try:
            self.molecules[name].update_conformation_map(conf_name, float(rmsd))
        except KeyError:
            warnings.warn(f"Attempt to update unknown molecule '{name}' with RMSD value '{rmsd}' "
                          f"for conformation '{conf_name}'")

    @staticmethod
    def _prepare_environment():
        error_messages: List[str] = []

        # input and template directories must exist as well as SB executable
        if not _Storage.input_dir.is_dir():
            error_messages.append(f"Input directory '{_Storage.input_dir}' was not found.")
        if not _Storage.templates_dir.is_dir():
            error_messages.append(f"Templates directory '{_Storage.templates_dir}' was not found.")
        if not _Storage.SB_executable_path.is_file():
            error_messages.append(f"SB executable '{_Storage.SB_executable_path}' not found.")
        if _Storage.crossplatform_runner and not which(_Storage.crossplatform_runner):
            error_messages.append(f"'{_Storage.crossplatform_runner}' command unavailable.")

        # output directory can already exist or be created
        try:
            _Storage.output_dir.mkdir(mode=0o777, parents=True, exist_ok=True)
        except (OSError, FileNotFoundError):
            error_messages.append("Output directory was not found and cannot be created.")

        if error_messages:
            raise ConfComparerError(" ".join(error_messages))


if __name__ == "__main__":
    # parse command line arguments
    parser = argparse.ArgumentParser(
        description="Determine conformation of carbon rings by comparing them to conformation "
                    "templates. Execution results will be stored in output dir, only brief "
                    "statistics are printed to standard output by default (this behaviour can "
                    "be altered by script arguments, see documentation of supported options).")
    parser.add_argument("-t", "--templates_dir", type=Path, default=_Storage.templates_dir,
                        help=f"path to directory containing templates (default: "
                             f"'{_Storage.templates_dir}')")
    parser.add_argument("-i", "--input_dir", type=Path, default=_Storage.input_dir,
                        help=f"path to directory containing input files (default: "
                             f"'{_Storage.input_dir}')")
    parser.add_argument("-o", "--output_dir", type=Path, default=_Storage.output_dir,
                        help=f"path to directory to be filled with results (will be created if "
                             f"does not exist, (default: '{_Storage.output_dir}')")
    parser.add_argument("-e", "--SB_executable_path", type=Path,
                        default=_Storage.SB_executable_path,
                        help=f"path to SiteBinder executable (default: "
                             f"'{_Storage.SB_executable_path}')")
    parser.add_argument("-d", "--delta_RMSD_tolerance", type=float, default=_Storage.tolerance,
                        help="highest acceptable RMSD deviation between conformation template and "
                             "tested molecule to be considered being in given conformation")
    parser.add_argument("-c", "--crossplatform_runner", type=str, default='',
                        help="command used to run SiteBinder executable (for Windows), typically "
                             "'mono' or 'wine' on Unix-like systems; if specified, analysis will "
                             "be executed as 'CROSS-PLATFORM_RUNNER PATH/TO/SB_EXECUTABLE "
                             "input.txt rmsd.csv pairings.csv'")
    parser.add_argument("-l", "--print_list", action="store_true",
                        help="print the execution results as list of molecules to standard output "
                             "(can be also found later in output directory)")
    parser.add_argument("-r", "--print_RMSD_chart", action="store_true",
                        help="print the execution results as RMSD chart in CSV format to standard "
                             "output (can be also found later in output directory)")
    parser.add_argument("-S", "--print_summary_off", action="store_true",
                        help="do not print the execution results as brief statistics to standard "
                             "output (still can be found later in output directory)")
    args = parser.parse_args()

    # load script parameters
    _Storage.templates_dir = args.templates_dir
    _Storage.input_dir = args.input_dir
    _Storage.output_dir = args.output_dir
    _Storage.SB_executable_path = args.SB_executable_path
    _Storage.crossplatform_runner = args.crossplatform_runner
    _Storage.tolerance = args.delta_RMSD_tolerance
    _Storage.print_list = args.print_list
    _Storage.print_RMSD_chart = args.print_RMSD_chart
    _Storage.print_summary = not args.print_summary_off

    analysis = Analysis()
    analysis.execute()
