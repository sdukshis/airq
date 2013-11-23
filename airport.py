"""
Airport class
"""

import re
import yaml
import math

from scipy import interpolate
import numpy
import pylab

from cs2cs import LatLonToUTMXY


class Airport(object):

    def __init__(self, data):
        self.data = data

    def get_sid(self, name):
        assert type(name) is str
        for sid in self.data['SID']:
            if sid['name'] == name:
                return Airline(direction='SID', data=sid)


class Airline(object):

    def __init__(self, direction, data):
        assert direction == 'SID' or direction == 'STAR'
        self.name = data['name']
        self.direction = direction
        self.RWY = data['RWY']
        self.inAngle = data['RWYAngle'] if direction == 'SID' else data['InAngle']
        self.outAngle = data['RWYAngle'] if direction == 'STAR' else data['OutAngle']
        self.RNAV = data['RNAV']

        self.points = self.convertRNAVtoUTMXY(self.RNAV)

    def convertRNAVtoUTMXY(self, RNAV):
        points = []
        for checkpoint in RNAV:
            wgs84string = checkpoint['wgs84']
            lat, lon = self.parse_wgs84(wgs84string)
            x, y = LatLonToUTMXY(lat, lon)
            points.append((x, y))

        return points

    def parse_wgs84(self, wgs84str):
        pattern = (r'^(?P<NS>N|S)(?P<LatA>\d{1,2})\s(?P<LatM>\d{1,2}\.\d)\s'
                   r'(?P<EW>E|W)(?P<LonA>\d{1,3})\s(?P<LonM>\d{1,2}\.\d)$')

        m = re.match(pattern, wgs84str)
        latsign = 1.0 if m.group('NS') == 'N' else -1.0
        lonsign = 1.0 if m.group('EW') == 'E' else -1.0
        lat = latsign*(int(m.group('LatA')) + float(m.group('LatM'))/60.0)
        lon = lonsign*(int(m.group('LonA')) + float(m.group('LonM'))/60.0)

        return (lat, lon)


if __name__ == "__main__":
    data = yaml.load(open('airports/AMS.air').read())
    AMS = Airport(data=data)
    AMGOD2X = AMS.get_sid(name='ANDIK 2E')
    # pprint(AMGOD2X.points)
    p = AMGOD2X.points

    def dist(p1, p2):
        return math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

    # for i in range(len(p) - 1):
    #     print(dist(p[i], p[i+1])/1000/1.852)

    x = map(lambda p: p[0], p)
    y = map(lambda p: p[1], p)

    tck,u = interpolate.splprep([x,y],s=0)
    unew = numpy.arange(0,1.01,0.01)
    out = interpolate.splev(unew,tck)

    
    pylab.figure()
    pylab.plot(x,y,'x',out[0],out[1])
    pylab.show()

    # print out
    # f = interpolate.interp1d(x, y)

    # xnew = np.arange(0,9, 0.1)
    # ynew = f(xnew)   # use interpolation function returned by `interp1d`
    # plt.plot(x, y, 'o', xnew, ynew, '-')
    # plt.show()




