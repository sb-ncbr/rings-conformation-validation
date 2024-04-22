from __future__ import annotations
from Molecules.Components.atom import Atom
from enum import Enum
from Utils.config import Config
from typing import Optional
from sys import stderr

from Utils.exceptions import InvalidSourceDataException


class MoleculeType(Enum):
    Undefined = 0
    Cyclopentane = 1
    Cyclohexane = 2
    Benzene = 3


class Conformation(Enum):
                    #    | CyH CyP Ben |
    Undefined = -1  #    |  x   x   x  |     Fallback option when no conformation was found
    Unanalysed = 0  #    |  x   x   x  | Default state of new molecule, no conformation found yet
                    #    |-------------|
    Flat = 1  #          |  x   x   x  |
    Half_Chair = 2  #    |  x   x   -  |
    Chair = 3  #         |  x   -   -  |
    Boat = 4  #          |  x   -   -  |
    Envelope = 5  #      |  -   x   -  |
    Tw_boat_right = 6  # |  x   -   -  |
    Tw_boat_left = 7  #  |  x   -   -  |


class Molecule:
    # A static variable containing data loaded from atom_names file on startup
    names: dict[str, list[set[str]]] = None
    # Config from which to load tolerances later
    config: Config = None
    # A list of all saved molecules
    molecules: list[Optional[Molecule]] = []

    @staticmethod
    def initialize(names: dict[str, list[set[str]]]):
        """
        Initializes the static parameters of the entire Molecule class
        on the first call.
        When processing data in multi-threaded mode, config gets set
        again within the worker class too.
        """
        if Molecule.names is None:
            Molecule.names = names

        if Molecule.config is None:
            Molecule.config = Config()

    def __init__(self, molecule_type: MoleculeType):
        # A list of all atoms within a given molecule
        self.atoms: list[Atom] = []

        # A list of possible conformation states the molecule can be in
        self.conformations: list[Conformation] = [Conformation.Unanalysed, Conformation.Undefined]

        self.conformation: Conformation = Conformation.Unanalysed
        self.molecule_type: MoleculeType = molecule_type
        self.is_valid: bool = True
        self.file_name: Optional[str] = None
        self.ligand: Optional[str] = None

        self.set_conformations()

    def print_statistics(self) -> None:
        """
        Calculates and prints out the summary of all the molecules that have
        been processed and their conformation has been decided.
        """
        print("SUMMARY\n-------")
        # Remove possible None-s from list from cases where file was ommited
        Molecule.molecules = [x for x in Molecule.molecules if x is not None]
        total = len(Molecule.molecules)
        for conf in self.conformations:
            count = sum([1 if x.conformation == conf else 0 for x in Molecule.molecules])
            percentage = (count / total) * 100
            percentage = int(percentage) if percentage % 1 == 0 else percentage
            print(f"{(conf.name.upper() + ':'):14}{count} ({percentage}%)")
        print(f"{'TOTAL:':14}{total}")

    def collect_data(self):
        dct = dict()
        Molecule.molecules = [x for x in Molecule.molecules if x is not None]
        # create empty dict with all possible conformations as keys
        for conf in self.conformations:
            dct[conf] = []

        # iterate over molecules, assign them to their conformation
        for molecule in Molecule.molecules:
            dct[molecule.conformation].append(molecule.file_name)

        rslt = dict()
        # rename names in dict so it's strings
        for conf in self.conformations:
            rslt[conf.name.lower()] = dct[conf]

        return rslt

    def __str__(self) -> str:
        return f"{self.file_name}: {self.conformation.name.upper().replace('_', ' ')}"

    def set_file_name(self, path: str) -> None:
        """
        Takes the path to the file and extracts the file name, setting appropriate
        field to its value
        """
        self.file_name = path.split("/")[-1].strip()

    def set_conformations(self) -> None:
        """
        Sets the correct order of possible conformations for each molecule
        which are later used for printing statistics.
        """
        match self.molecule_type:
            case MoleculeType.Cyclohexane:
                self.set_conf_cyclohexane()
            case MoleculeType.Cyclopentane:
                self.set_conf_cyclopentane()
            case MoleculeType.Benzene:
                self.set_conf_benzene()

    def set_conf_cyclohexane(self) -> None:
        self.conformations = ([Conformation.Boat, Conformation.Chair,
                              Conformation.Flat, Conformation.Half_Chair,
                              Conformation.Tw_boat_left, Conformation.Tw_boat_right]
                              + self.conformations)

    def set_conf_cyclopentane(self) -> None:
        self.conformations = [Conformation.Envelope, Conformation.Flat,
                              Conformation.Half_Chair] + self.conformations

    def set_conf_benzene(self) -> None:
        self.conformations = [Conformation.Flat] + self.conformations

    def get_atom_count(self):
        match self.molecule_type:
            case MoleculeType.Cyclohexane:
                return 6
            case MoleculeType.Benzene:
                return 6
            case MoleculeType.Cyclopentane:
                return 5

    def create_from_source(self, source: list[str]) -> None:
        """
        Creates atoms from the source file and detect ligand name.
        """
        for i in range(self.get_atom_count()):
            try:
                self.atoms.append(Atom(source[i + 1]))
            except InvalidSourceDataException:
                self.is_valid = False

        # detect and set residue name of the ligand if it exists
        self.is_valid = self.is_valid and self.atoms[0].residue_name is not None
        self.ligand = self.atoms[0].residue_name

    def validate_atoms(self) -> None:
        """
        Validates whether all the atoms of a molecule have their name
        present in the names file and then creates the ordering
        of atoms based on the order from the names file. In case of
        duplicate names within the same index of atom position, error
        out and cancel processing of this molecule as a result.
        """
        atom_count = self.get_atom_count()
        new_lst = [None for _ in range(atom_count)]
        for atom in self.atoms:
            for i in range(atom_count):
                if atom.name in self.names[self.ligand][i]:
                    if new_lst[i] is not None:
                        print(f"{atom.name} atom found twice!", file=stderr)
                        self.is_valid = False
                        return
                    new_lst[i] = atom
                    break
        if None in new_lst:
            print("Not all atoms were found")
            self.is_valid = False
            return
        self.atoms = new_lst
