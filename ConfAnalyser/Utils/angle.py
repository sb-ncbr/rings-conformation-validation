from ..Molecules.Components.geometries import Vector, Point
from math import degrees, atan2

def dihedral_angle(p1: Point, p2: Point, p3: Point, p4: Point) -> float:
    """
    Calculate a dihedral angle in degrees given four points in space.
    """
    u = Vector(p2, p1).cross(Vector(p2, p3))
    v = Vector(p3, p2).cross(Vector(p3, p4))
    unit_vector = Vector(p2, p3).normalize()

    return degrees(atan2(u.cross(v).dot(unit_vector), u.dot(v)))

