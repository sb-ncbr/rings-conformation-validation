from Molecules.Components.geometries import Plane
from Molecules.Components.rings import SixAtomRing
from Molecules.Components.molecule import MoleculeType, Conformation


class Benzene(SixAtomRing):
    def __init__(self, source_line: list[str], file_name: str):
        # Initialize the parent structure
        super().__init__(MoleculeType.Benzene)

        # Set the needed parameters
        self.source_file: list[str] = source_line
        self.set_file_name(file_name)

        try:
            self.create_from_source(source_line)
            self.validate_atoms()
            if self.is_valid:
                self.analyze()
        except Exception as e:  # should not happen but just in case, so we don't kill program
            print(e)

    def analyze(self):
        """
        Runs the analysis of the molecule, first finding the main plane of the
        molecule, after that it checks for the possible conformations to be true.
        """

        self.find_plane(self.config.b.t_flat_in)

        if self.has_plane and self.is_flat():
            self.conformation = Conformation.Flat
        else:
            self.conformation = Conformation.Undefined

    def is_flat(self):
        """
        Decides whether this molecule's conformation is flat.
        Flat conformation of molecule is determined by having all its atoms in one plane.
        """
        right_plane = Plane(self[0], self[1], self[3])
        left_plane = Plane(self[0], self[1], self[4])
        return (right_plane.is_on_plane(self[2], self.config.b.t_flat_in) and
                right_plane.is_on_plane(self[5], self.config.b.t_flat_in) and
                left_plane.is_on_plane(self[2], self.config.b.t_flat_in) and
                left_plane.is_on_plane(self[5], self.config.b.t_flat_in))
