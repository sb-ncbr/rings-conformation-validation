from enum import Enum, auto
from typing import List

class Ring(Enum):
    CYCLOHEXANE = auto()
    CYCLOPENTANE = auto()
    BENZENE = auto()
    OXANE=auto()
    OXOLANE=auto()

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
                return ["Twist"]
            


