#!/usr/bin/env python
from __future__ import print_function

from math import sin, cos, radians, pi


def main():
    theta = [ radians(x) for x in range(0, 360, 5) ]
    for t in theta:
        print('v %+0.2f %+0.2f +0.50' % ((10+7*sin(4*t))*cos(t), (10+7*cos(4*t))*sin(t)))
        print('vn 0 0 %+0.4f' % (t - pi,))

if __name__ == "__main__":
    main()
