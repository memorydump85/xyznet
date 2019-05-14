#!/usr/bin/env python
from __future__ import print_function

import fileinput


def main():
    infile = fileinput.input()
    line = infile.readline()
    assert line.rstrip('\n') == "OFF"

    num_vertices, num_faces, _ = [ int(x) for x in infile.readline().strip().split() ]
    for _ in range(num_vertices):
        print('v ' + infile.readline().strip())

    for _ in range(num_faces):
        fields = infile.readline().strip().split()
        num_sides = int(fields[0])
        v_indices = ( str(int(e)+1) for e in fields[1:num_sides+1] ) # OBJ indices are 1-based
        print('f ' + ' '.join(v_indices))

if __name__ == "__main__":
    main()
