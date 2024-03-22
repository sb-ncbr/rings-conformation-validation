from ..Molecules.Components.rings import SixAtomRing
from ..Molecules.Components.geometries import Plane
from ..Utils.angle import dihedral_angle
from ..Molecules.Components.molecule import MoleculeType, Conformation


class Cyclohexane(SixAtomRing):
    def __init__(self, source_line: list[str]):
        # Initialize the parent structure
        super().__init__(MoleculeType.Cyclohexane)

        # Set the needed parameters
        self.source_file: list[str] = source_line

        self.right_plane = None
        self.left_plane = None

        try:
            self.create_from_source(source_line)
            self.validate_atoms()
            if self.is_valid:
                self.analyze()
        except Exception as e:  # should not happen but just in case, so we don't kill program
            print(f"Error1: {e}")

    def analyze(self):
        """
        Runs the analysis of the molecule, first finding the main plane of the
        molecule, after that it checks for the possible conformations to be true.
        """

        self.find_plane(self.config.ch.t_in)

        self.right_plane = Plane(self[0], self[1], self[3])
        self.left_plane = Plane(self[0], self[1], self[4])

        if self.has_plane:
            if self.is_flat():
                self.conformation = Conformation.Flat
            elif self.is_half_chair():
                self.conformation = Conformation.Half_Chair
            elif self.is_chair():
                self.conformation = Conformation.Chair
            elif self.is_boat():
                 self.conformation = Conformation.Boat
            else:
                self.conformation = Conformation.Undefined
        else:
            if self.is_tw_boat_right():
                self.conformation = Conformation.Tw_boat_right
            elif self.is_tw_boat_left():
                self.conformation = Conformation.Tw_boat_left
            else:
                self.conformation = Conformation.Undefined

    def is_flat(self):
        """
        Decides whether this molecule's conformation is flat.
        Flat conformation of molecule is determined by having all its atoms in one plane.
        """
        return (self.right_plane.is_on_plane(self[2], self.config.ch.t_flat_in) and
                self.right_plane.is_on_plane(self[5], self.config.ch.t_flat_in) and
                self.left_plane.is_on_plane(self[2], self.config.ch.t_flat_in) and
                self.left_plane.is_on_plane(self[5], self.config.ch.t_flat_in))

    def is_half_chair(self) -> bool:
        """
        Decides whether this molecule's conformation is half chair.
        Half chair conformation of molecule is determined by having all but one atom in one plane.
        """
        right_dist = self.right_plane.true_distance_from(self[2])
        left_dist = self.right_plane.true_distance_from(self[5])
        return (self.right_plane.is_on_plane(self[2], self.config.ch.t_flat_in) !=
                self.right_plane.is_on_plane(self[5], self.config.ch.t_flat_in)) and \
                self.right_plane.is_on_plane(self[4], self.config.ch.t_flat_in) and \
                ((right_dist > self.config.ch.t_out) != (left_dist > self.config.ch.t_out))

    def is_chair(self) -> bool:
        """
        Decides whether the molecule's conformation is a chair.
        Chair conformation is determined by having two opposite atoms on
        opposite sides of the main molecule plane.
        """
        right_dist = self.right_plane.true_distance_from(self[2])
        left_dist = self.right_plane.true_distance_from(self[5])
        return right_dist > self.config.ch.t_out and \
               left_dist > self.config.ch.t_out and  \
               self.right_plane.are_opposite_side(self[2], self[5])

    def is_boat(self) -> bool:
        """
        Decides whether the molecule's conformation is a boat.
        Boat conformation is determined by having two opposite atoms on the
        same side of a plane.
        """
        right_dist = self.right_plane.true_distance_from(self[2])
        left_dist = self.right_plane.true_distance_from(self[5])
        return right_dist > self.config.ch.t_out and \
               left_dist > self.config.ch.t_out and \
               self.right_plane.are_same_side(self[2], self[5])


    def is_tw_boat_right(self) -> bool:
        """
        Determine whether the molecule's conformation is a twisted boat.
        Twisted boat conformation is determined by having no plane within the ring.

        The sight/left component is determined by having distance of atom 2
        from right plane be bigger than distance of atom 5 from left plane.
        """
        right_dist = self.right_plane.signed_distance_from(self[2])
        left_dist = self.left_plane.signed_distance_from(self[5])
        tw_angle = abs(dihedral_angle(self[1], self[3], self[4], self[0]))
        return ((tw_angle > self.config.ch.t_tw_boat_angle - self.config.ch.t_angle and
                tw_angle < self.config.ch.t_tw_boat_angle + self.config.ch.t_angle) and
                abs(right_dist) > self.config.ch.t_tw_out and
                abs(left_dist) > self.config.ch.t_tw_out and
                right_dist * left_dist > 0 and
                right_dist > left_dist)



    def is_tw_boat_left(self) -> bool:
        """
        Determine whether the molecule's conformation is a twisted boat.
        Twisted boat conformation is determined by having no plane within the ring.

        The sight/left component is determined by having distance of atom 2
        from right plane be smaller than distance of atom 5 from left plane.
        """
        right_dist = self.right_plane.signed_distance_from(self[2])
        left_dist = self.left_plane.signed_distance_from(self[5])
        tw_angle = abs(dihedral_angle(self[1], self[3], self[4], self[0]))

        return ((tw_angle > self.config.ch.t_tw_boat_angle - self.config.ch.t_angle and
                tw_angle < self.config.ch.t_tw_boat_angle + self.config.ch.t_angle) and
                abs(right_dist) > self.config.ch.t_tw_out and
                abs(left_dist) > self.config.ch.t_tw_out and
                right_dist * left_dist > 0 and
                right_dist < left_dist)
