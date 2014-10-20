"""
Airport class
"""

from __future__ import division, print_function
import re
import yaml
import math
from pprint import pprint

from scipy import interpolate
import numpy
# import pylab

from cs2cs import LatLonToUTMXY
from trajectory import trajectory
from qmodel import complexity2throughput

QUANTIM_STATE_LENGTH = 6790
SEPARATION_MINIMA = 9260

class Airport(object):

    def __init__(self, data):
        self.data = data

    def get_sid(self, name):
        assert type(name) is str
        for sid in self.data['SID']:
            if sid['name'] == name:
                return Airline(direction='SID', data=sid)

    def sids(self):
        for sid in self.data['SID']:
            yield Airline(direction='SID', data=sid)

    def stars(self):
        for star in self.data['STAR']:
            yield Airline(direction='STAR', data=star)

    def get_star(self, name):
        assert type(name) is str
        for star in self.data['STAR']:
            if star['name'] == name:
                return Airline(direction='STAR', data=star)


class Airline(object):

    def __init__(self, direction, data):
        assert direction == 'SID' or direction == 'STAR'
        self.name = data['name']
        self.direction = direction
        self.RWY = data['RWY']
        # self.inAngle = data['RWYAngle'] if direction == 'SID' else data['InAngle']
        # self.outAngle = data['RWYAngle'] if direction == 'STAR' else data['OutAngle']
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

    def interpolate(self):
        p = self.points
        v = self.get_velocity()
        return trajectory(p, v)

    def get_states(self, length):

        qstate_x = []
        qstate_y = []

        x, y = self.interpolate()

        x = map(lambda p: p[0], self.points)
        y = map(lambda p: p[1], self.points)

        curr_len = 0.0
        for i in range(len(x) - 1):
            curr_len += distance(x[i + 1], y[i + 1],
                                 x[i], y[i])

            if curr_len > length:
                qstate_x.append(x[i])
                qstate_y.append(y[i])
                curr_len = 0

        return (qstate_x, qstate_y)

    def get_points(self):
        x = map(lambda p: p[0], self.points)
        y = map(lambda p: p[1], self.points)

        tck, u = interpolate.splprep([x, y], s=0, k=2)
        unew = numpy.arange(0, 1.00, 0.01)
        out = interpolate.splev(unew, tck)

        return out

    def get_velocity(self):
        V = []
        for checkpoint in self.RNAV:
            Vkt = float(checkpoint['speed'])
            V.append(mph2mps(Vkt))

        return V

def mph2mps(v):
    """ Convert miles per hour to meters per sec """
    return (1.852 * v * 1000.0) / 3600.0

def distance(x1, y1, x2, y2):
    dx = x2 - x1
    dy = y2 - y1
    return math.sqrt(dx*dx + dy*dy)

def intersections(star, sid, separation_minima):
    star_x, star_y = star.get_states(QUANTIM_STATE_LENGTH)
    sid_x, sid_y = sid.get_states(QUANTIM_STATE_LENGTH)

    intersec = []

    for i in range(len(star_x)):
        for j in range(len(sid_x)):
            if distance(star_x[i], star_y[i],
                        sid_x[j], sid_y[j]) < separation_minima:
                intersec.append((i,j))

    return intersec


def topology_complexity(star, sid):
    p = []
    for i in intersections(star, sid, SEPARATION_MINIMA):
        l = len(star.get_states(QUANTIM_STATE_LENGTH)[0]) - i[0]+ 1
        k = i[1] + 1
        plk = min(l, k)/(l*k)
        p.append(plk)


    q = map(lambda i: 1 - i, p)

    kappa = 0.0
    for i in range(len(p)):
        prod = 1.0
        for j in range(len(q)):
            if j != i:
                prod *= q[j]
        kappa += p[i]*prod

    return kappa

def airport_throughput(airport):
    """ 
        Compute complexity and throughput for each pair
        of SID/STAR and average them to find out the
        airport throughput.
    """
    C = 0.0
    count = 0
    for sid in airport.sids():
        for star in airport.stars():
            q = topology_complexity(sid, star)
            Cr, Cm, Cp, P = complexity2throughput(1, q, 0.5)
            # print("%s vs %s: %g(%g)" % (sid.name, star.name, Cr/2, q))
            C += Cr/2
            count +=1

    print("total C = %g"  % (C/count))


if __name__ == "__main__":
    data = yaml.load(open('airports/SVO.air').read())
    SVO = Airport(data=data)

    airport_throughput(SVO)

    # pylab.figure()

    # for star in SVO.data['STAR']:
    #     p = SVO.get_star(star['name']).points
    #     v = SVO.get_star(star['name']).get_velocity()
    #     # pylab.plot(map(lambda p: p[0], p), map(lambda p: p[1], p), 'x')
    #     x, y = trajectory(p, v)
    #     pylab.plot(x, y)

    # for sid in SVO.data['SID']:
    #     p = SVO.get_sid(sid['name']).points
    #     v = SVO.get_sid(sid['name']).get_velocity()
    #     # pylab.plot(map(lambda p: p[0], p), map(lambda p: p[1], p), 'x')
    #     x, y = trajectory(p, v)
    #     pylab.plot(x, y)

    # pylab.show()
    # besta_01 = SVO.get_star('BESTA 01')
    # besta_4g = SVO.get_sid('BESTA 4G')

    # print("kappa = %g" % topology_complexity(besta_01, besta_4g))
    # pylab.figure()

    # points1 = besta_01.get_points()
    # points2 = besta_4g.get_points()

    # pylab.plot(points1[0], points1[1], 'r', points2[0], points2[1], 'g')
    
    # star_x, star_y = besta_01.get_states(QUANTIM_STATE_LENGTH)
    # sid_x, sid_y = besta_4g.get_states(QUANTIM_STATE_LENGTH)

    # pylab.plot(star_x, star_y, 'x')
    # pylab.plot(sid_x, sid_y, 'x')
    # for i in intersections(besta_01, besta_4g, SEPARATION_MINIMA):
    #     xs = [star_x[i[0]], sid_x[i[1]]]
    #     ys = [star_y[i[0]], sid_y[i[1]]]
    #     pylab.plot(xs, ys, 'b')
    #     print("%d - %d" % (i[0], i[1])) 
    # pylab.show()
    # for star in SVO.data['STAR']:
    #     p = SVO.get_star(star['name']).points

    #     x = map(lambda p: p[0], p)
    #     y = map(lambda p: p[1], p)

    #     tck, u = interpolate.splprep([x, y], s=0, k=2)
    #     unew = numpy.arange(0, 1.00, 0.01)
    #     out = interpolate.splev(unew, tck)

    #     qstate_x, qstate_y = SVO.get_star(star['name']).quantum_states(6790)

    #     pylab.plot(qstate_x, qstate_y, 'x', out[0], out[1])

    # for sid in SVO.data['SID']:
    #     p = SVO.get_sid(sid['name']).points

    #     x = map(lambda p: p[0], p)
    #     y = map(lambda p: p[1], p)

    #     tck, u = interpolate.splprep([x,y], s=0, k=2)
    #     unew = numpy.arange(0, 1.01, 0.01)
    #     out = interpolate.splev(unew, tck)

    #     qstate_x, qstate_y = SVO.get_sid(sid['name']).quantum_states(6790)
        
    #     pylab.plot(qstate_x, qstate_y, 'x', out[0], out[1])

    # pylab.show()
    # AMGOD2X = AMS.get_sid(name='ANDIK 2E')
    # # pprint(AMGOD2X.points)
    # p = AMGOD2X.points

    # def dist(p1, p2):
    #     return math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

    # # for i in range(len(p) - 1):
    # #     print(dist(p[i], p[i+1])/1000/1.852)

    # x = map(lambda p: p[0], p)
    # y = map(lambda p: p[1], p)

    # tck,u = interpolate.splprep([x,y],s=0)
    # unew = numpy.arange(0,1.01,0.01)
    # out = interpolate.splev(unew,tck)

    
    # pylab.figure()
    # pylab.plot(x,y,'x',out[0],out[1])
    # pylab.show()

    # print out
    # f = interpolate.interp1d(x, y)

    # xnew = numpy.arange(0,9, 0.1)
    # ynew = f(xnew)   # use interpolation function returned by `interp1d`
    # plt.plot(x, y, 'o', xnew, ynew, '-')
    # plt.show()
