#!/bin/env python

import sys

from array import array

class Point:
    def __init__(self, x, y, z, nx, ny, nz):
        self.x = x
        self.y = y
        self.z = z
        self.nx = nx
        self.ny = ny
        self.nz = nz

    def print_point(self):
        print('(%(x)f, %(y)f, %(z)f) (%(nx)f, %(ny)f, %(nz)f)' % self)

def string_to_point(string):
    values = [float(x) for x in string.split()]
    return Point(values[0], values[1], values[2], values[3], values[4], values[5])

def parse_text(filename):
    f = open(filename, 'r')
    points = [string_to_point(line) for line in f]
    f.close()
    return points

def print_binary(filename, points):
    f = open(filename, 'wb')
    float_array = array('f')
    for p in points:
        float_array.fromlist([p.x, p.y, p.z, p.nx, p.ny, p.nz])
    float_array.tofile(f)
    f.close()

if len(sys.argv) != 3:
    print('USAGE: ' + sys.argv[0] + ' [input-file] [output-file]')
    exit(1)

file_from = sys.argv[1]
file_to = sys.argv[2]

print('Reading from ' + file_from + ' and writing to ' + file_to)

points = parse_text(file_from)
print_binary(file_to, points)
