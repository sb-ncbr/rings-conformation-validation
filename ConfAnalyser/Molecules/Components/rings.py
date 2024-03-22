from ...Molecules.Components.molecule import Molecule, MoleculeType, Conformation
from ...Molecules.Components.geometries import Plane


class Ring(Molecule):
    def __init__(self, molecule_type: MoleculeType):
        super().__init__(molecule_type)
        self.conformation = Conformation.Undefined
        self.has_plane = False
        self.begin = 0

    def __getitem__(self, item):
        """
        A getter so we can in code get atoms on specific index by simply
        using self[index]. This way we can avoid having to always include
        calculation of actual index in regard to `begin` field. This will
        calculate this for us every time so the code can be cleaner
        """
        if isinstance(item, int):
            return self.atoms[(self.begin + item) % self.get_atom_count()]
        return None


class SixAtomRing(Ring):
    def __init__(self, molecule_type: MoleculeType):
        super().__init__(molecule_type)

    def find_plane(self, tolerance: float, dist1: int = 1, dist2: int = 3, dist3: int = 4) -> bool:
        """
        Find the best suitable plane to start working with and set the `begin` parameter
        of the molecule to be used in other conformation validators.
        Also decide whether the molecule even has any valid plane at all.
        """
        distance = float("inf")
        self.has_plane = False

        for i in range(6):
            plane = Plane(self.atoms[i % 6], self.atoms[(i + dist1) % 6], self.atoms[(i + dist2) % 6])
            dist_from_plane = abs(plane.true_distance_from(self.atoms[(i + dist3) % 6]))
            if dist_from_plane < distance:
                self.begin = i
                distance = dist_from_plane

                if plane.is_on_plane(self.atoms[(i + dist3) % 6], tolerance):
                    self.has_plane = True

        return self.has_plane


class FiveAtomRing(Ring):
    def __init__(self, molecule_type: MoleculeType):
        super().__init__(molecule_type)

    def find_plane(self, tolerance: float, dist1: int = 1, dist2: int = 2, dist3: int = 3) -> bool:
        """
        Find the best suitable plane to start working with and set the `begin` parameter
        of the molecule to be used in other conformation validators.
        Also decide whether the molecule even has any valid plane at all.
        """
        distance = float("inf")

        for i in range(5):
            plane = Plane(self.atoms[i % 5], self.atoms[(i + dist1) % 5], self.atoms[(i + dist2) % 5])
            dist_from_plane = abs(plane.true_distance_from(self.atoms[(i + dist3) % 5]))
            if dist_from_plane < distance:
                self.begin = i
                distance = dist_from_plane

                if plane.is_on_plane(self.atoms[(i + dist3) % 5], tolerance):
                    self.has_plane = True

        return self.has_plane
