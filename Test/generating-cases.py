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
    def mult(self, val):
        return Vector(self.x * val, self.y * val, self.z * val)
    def sum(self, vec):
        return Vector(self.x + vec.x, self.y + vec.y, self.z + vec.z)
    def normalize(self):
        return self.mult(1 / self.norm())
    def random_pos(self, eps):
        vec = random_vector().mult(eps)
        return self.sum(vec)
    def to_radial(self):
        r = self.norm()
        theta = math.acos(self.z / r)
        phi = math.atan2(self.y, self.x)
        return (r, theta, phi)
    def random_orient(self, eps):
        (r, theta, phi) = self.to_radial()
        (newtheta, newphi) = random_radials(theta, phi, eps)
        return from_radial(r, newtheta, newphi)

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

def points_to_text(points):
    for p in points:
        p.print()

def generate_sphere(counts, radius, pos_eps, angle_eps):
    points = []
    for i in range(0, counts):
        norm = random_vector()
        pos = norm.mult(radius)
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

def do_sphere():
    points_to_text(generate_sphere(1000, 100, 5, math.pi / 10))

def do_cube():
    points_to_text(generate_cube(3000, 200, 5, math.pi / 10))

if len(sys.argv) != 2:
    print('USAGE: ', sys.argv[0], ' [cube|sphere]')
    sys.exit(1)

if sys.argv[1] == 'sphere':
    do_sphere()
elif sys.argv[1] == 'cube':
    do_cube()
