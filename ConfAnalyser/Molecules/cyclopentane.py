from Molecules.Components.geometries import Plane
from Molecules.Components.rings import FiveAtomRing
from Molecules.Components.molecule import MoleculeType, Conformation
from Utils.angle import dihedral_angle


class Cyclopentane(FiveAtomRing):
    def __init__(self, source_line: list[str]):
        # Initialize the parent structure
        super().__init__(MoleculeType.Cyclopentane)

        # Set the needed parameters
        self.source_file: list[str] = source_line

        self.basePlane = None

        try:
            self.create_from_source(source_line)
            self.validate_atoms()
            if self.is_valid:
                self.analyze()
        except Exception as e:  # should not happen but just in case, so we don't kill program
            print(f"Error 2: {e}")

    def analyze(self):
        """
        Runs the analysis of the molecule, first finding the main plane of the
        molecule, after that it checks for the possible conformations to be true.
        """

        self.find_plane(self.config.cp.t_in)

        self.basePlane = Plane(self[0], self[1], self[2])

        if self.has_plane:
            if self.is_flat():
                self.conformation = Conformation.Flat
            elif self.is_envelope():
                self.conformation = Conformation.Envelope
            else:
                self.conformation = Conformation.Undefined
        else:
            if self.is_half_chair():
                self.conformation = Conformation.Half_Chair
            else:
                self.conformation = Conformation.Undefined


    def is_flat(self) -> bool:
        """
        Decides whether this molecule's conformation is flat.
        Flat conformation of molecule is determined by having all its atoms in one plane.
        """
        return (self.basePlane.is_on_plane(self[3], self.config.cp.t_in) and
                self.basePlane.is_on_plane(self[4], self.config.cp.t_in))

    def is_envelope(self) -> bool:
        """
        Decides whether this molecule's conformation is envelope.
        Envelope conformation of molecule is determined by having all but one atom on a same plane.
        Checked that with adding condition for if atom 3 is on plane and results
        did not change at all, probably extremely improbable
        """
        return self.basePlane.true_distance_from(self[4]) > self.config.cp.t_out

    def is_half_chair(self):
        """
        Decides whether this molecule's conformation is half chair.
        Twist conformation of molecule is determined by having no plane within the ring.
        """
        left_plane = Plane(self[0], self[1], self[3])
        right_plane = Plane(self[0], self[2], self[3])
        left_distance = left_plane.signed_distance_from(self[4])
        right_distance = right_plane.signed_distance_from(self[4])
        twist_angle = dihedral_angle(self[0], self[1], self[2], self[3])

        return ((abs(abs(twist_angle) - self.config.cp.t_tw_boat_angle) < self.config.cp.t_angle)
                and (abs(right_distance) > self.config.cp.t_tw_out)
                and (abs(left_distance) > self.config.cp.t_tw_out)
                and (right_distance * left_distance > 0))
