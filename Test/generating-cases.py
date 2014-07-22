#!/bin/env python

import math
import random
import sys

# eps in radians
def random_radials(theta, phi, eps):
    newphi = phi + (random.random() - 0.5) * 2 * eps
    newtheta = theta + (random.random() - 0.5) * 2 * eps
    if newtheta < 0:
        newtheta = -newtheta
        newphi = newphi + math.pi
    while newtheta > 2 * math.pi:
        newtheta = newtheta - 2 * math.pi
    if newtheta > math.pi:
        newtheta = 2 * math.pi - newtheta
        newphi = newphi + math.pi
    while newphi > 2 * math.pi:
        newphi = newphi - 2 * math.pi
    return (newtheta, newphi)

class Vector:
    def __init__(self, x = 0, y = 0, z = 0):
        self.x = x
        self.y = y
        self.z = z
    def norm2(self):
        return self.x * self.x + self.y * self.y + self.z * self.z
    def norm(self):
        return math.sqrt(self.norm2())
    def __mul__(self, val):
        return Vector(self.x * val, self.y * val, self.z * val)
    def __add__(self, vec):
        return Vector(self.x + vec.x, self.y + vec.y, self.z + vec.z)
    def __sub__(self, vec):
        return self + (-vec)
    def __neg__(self):
        return Vector(-self.x, -self.y, -self.z)
    def normalize(self):
        return self * (1 / self.norm())
    def random_pos(self, eps):
        vec = random_vector() * eps
        return self + vec
    def to_radial(self):
        r = self.norm()
        theta = math.acos(self.z / r)
        phi = math.atan2(self.y, self.x)
        return (r, theta, phi)
    def random_orient(self, eps):
        (r, theta, phi) = self.to_radial()
        (newtheta, newphi) = random_radials(theta, phi, eps)
        return from_radial(r, newtheta, newphi)
    def rotateX(self, angle):
        return Vector(self.x,
                self.y * math.cos(angle) + self.z * math.sin(angle),
                -self.y * math.sin(angle) + self.z * math.cos(angle))

def from_radial(r, theta, phi):
    z = r * math.cos(theta)
    x = r * math.sin(theta) * math.cos(phi)
    y = r * math.sin(theta) * math.sin(phi)
    return Vector(x, y, z)

def random_vector():
    return Vector(random.random() - 0.5, random.random() - 0.5, random.random() - 0.5).normalize()

class Point:
    def __init__(self, pos = Vector(), norm = Vector()):
        self.pos = pos
        self.norm = norm
    def print(self):
        print(self.pos.x, self.pos.y, self.pos.z, self.norm.x, self.norm.y, self.norm.z)
    def rotateX(self, angle):
        return Point(self.pos.rotateX(angle), self.norm.rotateX(angle))
    def shift(self, v):
        return Point(self.pos + v, self.norm)

def points_to_text(points):
    for p in points:
        p.print()

def generate_sphere(counts, radius, pos_eps, angle_eps):
    points = []
    for i in range(0, counts):
        norm = random_vector()
        pos = norm * radius
        points.append(Point(pos.random_pos(pos_eps), norm.random_orient(angle_eps)))
    return points

def generate_cube(counts, size, pos_eps, angle_eps):
    points = []
    for i in range(0, counts):
        px = (random.random() - 0.5) * size
        py = (random.random() - 0.5) * size
        pt = Point()
        if i % 6 == 0:
            pt = Point(Vector(px, py, size / 2), Vector(0, 0, 1))
        elif i % 6 == 1:
            pt = Point(Vector(px, py, -size / 2), Vector(0, 0, -1))
        elif i % 6 == 2:
            pt = Point(Vector(px, size / 2, py), Vector(0, 1, 0))
        elif i % 6 == 3:
            pt = Point(Vector(px, -size / 2, py), Vector(0, -1, 0))
        elif i % 6 == 4:
            pt = Point(Vector(size / 2, px, py), Vector(1, 0, 0))
        elif i % 6 == 5:
            pt = Point(Vector(-size / 2, px, py), Vector(-1, 0, 0))
        points.append(Point(pt.pos.random_pos(pos_eps), pt.norm.random_orient(angle_eps)))
    return points

def circle_segment_radius(r, phispan, phi):
    return r / (math.cos(phi) + math.sin(phi) * math.tan(phispan / 2))

def generate_cylinder(count, cap_count, zstart, zend, radius, phistart, phiend,
        pos_eps, pos_sigma, angle_eps):
    points = []
    for i in range(0, count):
        phi = random.random() * (phiend - phistart) + phistart
        norm = Vector(math.cos(phi), math.sin(phi), 0)
        pos = norm * (radius + random.gauss(0, pos_sigma) * pos_eps)
        pos.z = random.random() * (zend - zstart) + zstart
        points.append(Point(pos, norm.random_orient(angle_eps)))
    for i in range(0, 2 * cap_count):
        phi = random.random() * (phiend - phistart) + phistart
        b = circle_segment_radius(radius, phiend - phistart, phi - phistart)
        r = random.random() * (radius - b) + b
        theta = math.pi / 2
        pos = from_radial(r, theta, phi)
        if i < cap_count:
            pos.z = zstart
            norm = Vector(0, 0, -1)
        else:
            pos.z = zend
            norm = Vector(0, 0, 1)
        points.append(Point(pos, norm))
    for i in range(0, cap_count):
        phi = random.random() * (phiend - phistart) + phistart
        b = circle_segment_radius(radius, phiend - phistart, phi - phistart)
        theta = math.pi / 2
        pos = from_radial(b, theta, phi)
        pos.z = random.random() * (zend - zstart) + zstart
        norm = from_radial(1, theta, (phistart + phiend) / 2 + math.pi)
        points.append(Point(pos, norm))
    return points

def generate_cylinders(count, cap_count, length, radius, phistart, phispan, angle,
        pos_eps, pos_sigma, angle_eps):
    zstart = radius / math.tan(angle / 2)
    zend = zstart + length
    phiend = phistart + phispan
    cyl1 = generate_cylinder(int(count / 2), cap_count, zstart, zend, radius,
            phistart, phiend, pos_eps, pos_sigma, angle_eps)
    cyl2 = generate_cylinder(int(count / 2), cap_count, zstart, zend, radius,
            phistart + math.pi, phiend + math.pi, pos_eps, pos_sigma, angle_eps)
    cyls = []
    shiftV = Vector(-radius / 2, radius / 2, -zstart + length / 2)
    for p in cyl1:
        cyls.append(p.shift(shiftV))
    for p in cyl2:
        cyls.append(p.rotateX(angle).shift(shiftV))
    return cyls

def do_sphere():
    points_to_text(generate_sphere(1000, 100, 5, math.pi / 10))

def do_cube():
    points_to_text(generate_cube(3000, 200, 5, math.pi / 10))

def do_cylinders():
    points_to_text(generate_cylinders(30000, 10000, 50, 100, 3 * math.pi / 8, math.pi / 4, math.pi / 18,
        6, 0.2, math.pi / 100))

if len(sys.argv) != 2:
    print('USAGE: ', sys.argv[0], ' [cube|sphere|cylinders]')
    sys.exit(1)

if sys.argv[1] == 'sphere':
    do_sphere()
elif sys.argv[1] == 'cube':
    do_cube()
elif sys.argv[1] == 'cylinders':
    do_cylinders()
