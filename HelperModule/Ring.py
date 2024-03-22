from enum import Enum, auto

class Ring(Enum):
    CYCLOHEXANE = auto()
    CYCLOPENTANE = auto()
    BENZENE = auto()

    @property
    def atom_number(self) -> int:
        match self:
            case Ring.CYCLOHEXANE:
                return 6
            case Ring.BENZENE:
                return 6
            case Ring.CYCLOPENTANE:
                return 5

    @property
    def name_substring(self) -> str | None:
        match self:
            case Ring.CYCLOHEXANE:
                return 'cyclohex'
            case Ring.BENZENE:
                return 'benz;phen'
            case Ring.CYCLOPENTANE:
                return 'cyclopent'
            case _:
                return None

    @property
    def pattern_query(self) -> str | None:
        match self:
            case Ring.CYCLOHEXANE:
                return "Rings(6 * ['C'])"
            case Ring.BENZENE:
                return "Rings(6 * ['C'])"
            case Ring.CYCLOPENTANE:
                return "Rings(5 * ['C'])"
            case _:
                return None
