from __future__ import print_function, division

import argparse
import re
from math import sin, cos, radians

from cs2cs import LatLonToUTMXY

def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("angle", type=int)
    argparser.add_argument("radius", type=float)
    argparser.add_argument("wgs84")

    args = argparser.parse_args()

    lat, lon = parse_wgs84(args.wgs84)

    x, y = LatLonToUTMXY(lat, lon)

    xp, yp = intersect(x, y, miles2km(args.radius), args.angle)

    print(xp, yp)

def parse_wgs84(wgs84str):
    pattern = (r'^(?P<NS>N|S)(?P<LatA>\d{1,2})\s(?P<LatM>\d{1,2}\.\d)\s'
               r'(?P<EW>E|W)(?P<LonA>\d{1,3})\s(?P<LonM>\d{1,2}\.\d)$')

    m = re.match(pattern, wgs84str)
    latsign = 1.0 if m.group('NS') == 'N' else -1.0
    lonsign = 1.0 if m.group('EW') == 'E' else -1.0
    lat = latsign*(int(m.group('LatA')) + float(m.group('LatM'))/60.0)
    lon = lonsign*(int(m.group('LonA')) + float(m.group('LonM'))/60.0)

    return (lat, lon)


def intersect(x, y, radius, angle):
    return (x + radius*cos(radians(angle)),
            y + radius*sin(radians(angle)))


def miles2km(m):
    return m/1.852
if __name__ == "__main__":
    main()

