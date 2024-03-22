from __future__ import annotations
from math import sqrt
from typing import Union


class Point:
    def __init__(self, x: float = 0, y: float = 0, z: float = 0):
        """
        A class representing a single point in a 3-dimensional space, consisting of
        three float numbers.
        """
        self.x = x
        self.y = y
        self.z = z

    def distance_from(self, other: Point) -> float:
        """
        Calculates the distance of this point from the other point
        Formula: sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
        """
        return abs(sqrt((other.x - self.x) ** 2 + (other.y - self.y) ** 2 + (other.z - self.z) ** 2))

    def __str__(self):
        return f"Point: <{self.__repr__()}>"

    def __repr__(self):
        return f"{round(self.x, 3)}, {round(self.y, 3)}, {round(self.z, 3)}"


class Vector(Point):
    def __init__(self, x: Union[float, Point] = 0.0, y: Union[float, Point] = 0.0, z: float = 0.0):
        # Creating vector from two points
        if isinstance(x, Point) and isinstance(y, Point) and isinstance(z, float):
            p1, p2 = x, y
            super().__init__(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z)

        # Creating vector from a single point
        elif isinstance(y, float) and isinstance(x, Point) and isinstance(z, float):
            p1 = x
            super().__init__(p1.x, p1.y, p1.z)

        # Creatubg vector from three coordinates
        elif isinstance(x, float) and isinstance(y, float) and isinstance(z, float):
            super().__init__(x, y, z)

    def from_point(self, point: Point) -> Vector:
        """
        Create a vector from a single point in space
        """
        self.x = point.x
        self.y = point.y
        self.z = point.z

        return self

    def from_points(self, point1: Point, point2: Point) -> Vector:
        """
        Create a vector between the two provided points
        """
        self.x = point1.x - point2.x
        self.y = point1.y - point2.y
        self.z = point1.z - point2.z

        return self

    def __add__(self, other: Vector) -> Vector:
        return Vector(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other: Vector) -> Vector:
        return Vector(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, other: float) -> Vector:
        return Vector(self.x * other, self.y * other, self.z * other)

    def length(self) -> float:
        """
        Calculates the length/size of this vector
        """
        return sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

    def normalize(self) -> Vector:
        """
        Calculates a normalized, unit vector of this vector and returns it
        """
        length = self.length()
        return Vector(self.x / length, self.y / length, self.z / length)

    def dot(self, other: Vector) -> float:
        """
        Calculates a dot product between this and the other vector
        """
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cross(self, other: Vector) -> Vector:
        """
        Calculate a cross product between this and the other vector
        """
        return Vector(self.y * other.z - self.z * other.y,
                      self.z * other.x - self.x * other.z,
                      self.x * other.y - self.y * other.x)

    def __str__(self):
        return f"Vector ({self.__repr__()})"


class Plane:
    def __init__(self, point1: Point, point2: Point, point3: Point):
        self.u = Vector(point3, point2)
        self.v = Vector(point2, point1)

        self.normal = self.u.cross(self.v)
        self.d = -(self.normal.x * point1.x + self.normal.y * point1.y + self.normal.z * point1.z)
        self.a = self.normal.x
        self.b = self.normal.y
        self.c = self.normal.z

    def true_distance_from(self, point: Point) -> float:
        """
        Calculates a distance of a given point from this plane.
        """
        return abs(self.signed_distance_from(point))

    def signed_distance_from(self, point: Point) -> float:
        """
        Calculates a distance of a given point from this plane.

        The usual equation for distance between point and a plane actually uses
        absolute value inside of it to remove negative values as those make no
        sense within what we consider to be a distance. This function
        calculates a "partial" distance, as in it does not use absolute value
        and the result thus can be negative.

        Thanks to this we can determine on which side of the plane the point is
        located as all the points on one side of plane will have same sign and
        all the points on the other side of the plane will have opposite sign.
        """
        size = sqrt(self.a ** 2 + self.b ** 2 + self.c ** 2)
        if size == 0:  # shouldn't happen if data is correct
            return 0
        return ((self.a * point.x + self.b * point.y + self.c * point.z + self.d) / size)

    def is_on_plane(self, point: Point, tolerance: float) -> bool:
        """
        Decides whether a given point is located on this plane.
        Threshold determined by passed `tolerance` parameter
        """
        return self.true_distance_from(point) <= tolerance

    def are_opposite_side(self, point1: Point, point2: Point) -> bool:
        """
        Decides whether two provided points are on opposite sides of this plane.

        Calculate signed distances of both points from this plane, multiply
        then those distance. Since they were signed, if the result is
        negative, that means one of the distances was negative and
        other was positive which then means that the points are located on the
        opposite sides of this plane.
        """
        distance_1 = self.signed_distance_from(point1)
        distance_2 = self.signed_distance_from(point2)
        return distance_1 * distance_2 < 0

    def are_same_side(self, point1: Point, point2: Point) -> bool:
        """
        Decides whether two provided points are on the same side of this plane.
        """
        return not self.are_opposite_side(point1, point2)

    def __str__(self):
        return f"Plane of a, b, c, d = {self.a, self.b, self.c, self.d}"
