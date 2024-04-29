from Utils.exceptions import *
from Molecules.Components.geometries import Point


class Atom(Point):
    def __init__(self, source: str) -> None:
        """
        Constructor for creating the atom structure from a line of text
        loaded from the source .pdb file.

        :param source: single line representing atom from which object is created

        Source: https://www.biostat.jhsph.edu/~iruczins/teaching/260.655/links/pdbformat.pdf

        The structure follows pattern noted in source above ^:
        COLUMN  ALIGN   TYPE    DATA
        1-4     left    str     Either "ATOM" or "HETATM"
        7-11    right   int     Atom serial number
        13-16   left    str     Atom name
        17      ----    char    Alternate location indicator
        18-20   right   str     Residue name
        22      ----    char    Chain identifier
        23-26   right   int     Residue sequence number
        27      ----    char    Code for insertions of residues
        31-38   right   float   X orthogonal Angstrom coordinate
        39-46   right   float   Y orthogonal Angstrom coordinate
        47-54   right   float   Z orthogonal Angstrom coordinate
        55-60   right   float   Occupancy
        61-66   right   float   Temperature factor
        73-76   left    str     (optional) Segment identifier
        77-78   right   str     Element symbol
        79-80   ----    str     (optional) Charge
        """
        super().__init__()
        self.atom_type = None               # 1-4
        self.serial_number = None           # 7-11
        self.name = None                    # 13-16
        self.alt_loc_ident = None           # 17
        self.residue_name = None            # 18-20
        self.chain_identifier = None        # 22
        self.residue_sequence_num = None    # 23-26
        self.res_insert_code = None         # 27
        self.x = None                       # 31-38
        self.y = None                       # 39-46
        self.z = None                       # 47-54
        self.occupancy = None               # 55-60
        self.temperature_factor = None      # 61-66
        self.segment_identifier = None      # 73-76
        self.element_symbol = None          # 77-78
        self.charge = None                  # 79-80

        self.is_valid = True
        self.source = source

        try:
            self.create_from_source(source)
        except InvalidSourceDataException as e:
            self.is_valid = False
            print(f"[ERROR]: Atom - invalid source data: {e}")

    def create_from_source(self, source: str) -> None:
        """
        Parse a single line from the source .pdb file, parsing information
        from it and validating data validity, setting all the fields of this
        class object.

        Throw `InvalidSourceDataException` in case of an error.

        :param source: a line containing data from the source .pdb file
                       representing atom in the molecule structure
        :return: None
        """
        atom_type = source[0:6]
        if atom_type != "ATOM" and atom_type != "HETATM":
            raise InvalidSourceDataException(f"'{atom_type}' is not a valid atom type")
        self.atom_type = atom_type

        atom_serial_name = source[6:11].lstrip()
        if not atom_serial_name.isnumeric():
            raise InvalidSourceDataException(f"'{atom_serial_name}' is not a valid atom serial number")
        self.serial_number = int(atom_serial_name)

        atom_name = source[12:16].strip()
        if not atom_name:
            raise InvalidSourceDataException("Atom name is empty")
        self.name = atom_name

        self.alt_loc_ident = source[16] if source[16].strip() else None

        self.residue_name = source[17:20].strip()

        self.chain_identifier = source[21]

        residue_sequence_num = source[22:26].lstrip()
        try:
            self.residue_sequence_num = int(residue_sequence_num)
        except ValueError:
            raise InvalidSourceDataException(f"'{residue_sequence_num}' is not a valid residue sequence number")

        self.res_insert_code = source[26]

        coordinate_X = source[30:38].strip()
        try:
            self.x = float(coordinate_X)
        except ValueError:
            raise InvalidSourceDataException(f"'{coordinate_X}' is not a valid float number")

        coordinate_Y = source[38:46].strip()
        try:
            self.y = float(coordinate_Y)
        except ValueError:
            raise InvalidSourceDataException(f"'{coordinate_Y}' is not a valid float number")

        coordinate_Z = source[46:54].strip()
        try:
            self.z = float(coordinate_Z)
        except ValueError:
            raise InvalidSourceDataException(f"'{coordinate_Z}' is not a valid float number")

        occupancy = source[54:60].strip()
        try:
            self.occupancy = float(occupancy)
        except ValueError:
            raise InvalidSourceDataException(f"'{occupancy}' is not a valid float number")

        temperature_factor = source[60:66].strip()
        try:
            self.temperature_factor = float(temperature_factor)
        except ValueError:
            raise InvalidSourceDataException(f"'{temperature_factor}' is not a valid float number")

        segment_identifier = source[72:76].strip()
        self.segment_identifier = segment_identifier if segment_identifier else None

        element_symbol = source[76:78].strip()
        self.element_symbol = element_symbol if element_symbol else None

        charge = source[78:80].strip()
        self.charge = charge if charge else None

    def __str__(self):
        output = (f"Atom type: '{self.atom_type}'\n"
                  f"Atom serial number: {self.serial_number}\n"
                  f"Atom name: '{self.name}'\n"
                  f"Alternate location indicator: '{self.alt_loc_ident}'\n"
                  f"Residue name: '{self.residue_name}'\n"
                  f"Chain identifier: '{self.chain_identifier}'\n"
                  f"Residue sequence number: {self.residue_sequence_num}\n"
                  f"Code for insertions of residues: '{self.res_insert_code}'\n"
                  f"X: {self.x}\n"
                  f"Y: {self.y}\n"
                  f"Z: {self.z}\n"
                  f"Occupancy: {self.occupancy}\n"
                  f"Temperature factor: {self.temperature_factor}\n"
                  f"Segment identifier: '{self.segment_identifier}'\n"
                  f"Element symbol: '{self.element_symbol}'\n"
                  f"Charge: '{self.charge}'")
        output += f"Atom {self.name} - <{self.x}, {self.y}, {self.z}>"
        return output

