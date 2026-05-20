from enum import Enum, auto
from typing import List

class Ring(Enum):
    CYCLOHEXANE = auto()
    CYCLOPENTANE = auto()
    BENZENE = auto()
    OXANE=auto()
    OXOLANE=auto()
    PYRROLIDINE = auto()    # all single bonds
    PYRROLE = auto()    # aromatic
    PYRROLINE = auto()  # one double bond

    @property
    def atoms(self) -> str:
        match self:
            case Ring.CYCLOHEXANE:
                return "C*6"
            case Ring.BENZENE:
                return "C*6"
            case Ring.CYCLOPENTANE:
                return "C*5"
            case Ring.OXANE:
                return "C*5-O*1"
            case Ring.OXOLANE:
                return "C*4-O*1"
            case Ring.PYRROLIDINE:
                return "C*4-N*1"
            case Ring.PYRROLE:
                return "C*4-N*1"
            case Ring.PYRROLINE:
                return "C*4-N*1"

    @property
    def atom_number(self) -> int:
        match self:
            case Ring.CYCLOHEXANE:
                return 6
            case Ring.BENZENE:
                return 6
            case Ring.CYCLOPENTANE:
                return 5
            case Ring.OXANE:
                return 6
            case Ring.OXOLANE:
                return 5
            case Ring.PYRROLIDINE:
                return 5
            case Ring.PYRROLE:
                return 5
            case Ring.PYRROLINE:
                return 5
            
    @property
    def favourable_confs(self) -> List[str]:
        match self:
            case Ring.CYCLOHEXANE:
                return ["Chair"]
            case Ring.BENZENE:
                return ["Flat"]
            case Ring.CYCLOPENTANE:
                return ["Envelope", "Half chair"]
            case Ring.OXANE:
                return ["Chair"]
            case Ring.OXOLANE:
                return ["Twist", "Envelope"]
            


