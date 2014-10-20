# -*- coding: utf-8 -*-

from math import fabs, sin, cos, tan, atan2, sqrt, pi
from numpy import sign, arange

# import pylab

g = 9.80665

def trajectory(p, V):
    assert type(p) is list
    assert len(p) >= 3
    assert type(V) is list
    x = [p[0][0]]
    y = [p[0][1]]


    for i in range(len(p) - 2):
        # FIXME: продумать угол тангажа
        R = LUR_R(V[i + 1], pi/12);
        xt, yt = plot(p[i], p[i + 1], p[i + 2], R)
        # x0 = x[-1]
        # y0 = y[-1]
        # x1 = xt[0]
        # y1 = yt[0]
        # for t in arange(0, 1, 0.001):
        #     x.append(x0*(1 - t) + x1*t)
        #     y.append(y0*(1 - t) + y1*t)
        x.extend(xt)
        y.extend(yt)
        
    #Last segment to RunWay
    x0 = x[-1]
    y0 = y[-1]
    x1 = p[-1][0]
    y1 = p[-1][1]
    for t in arange(0, 1, 0.001):
        x.append(x0*(1 - t) + x1*t)
        y.append(y0*(1 - t) + y1*t)       
    return x, y

def plot(p1, p2, p3, R):
    b, e, c = LUR(p1, p2, p3, R)

    psi = angle_line(c, b)
    phi = angle_line(c, e)

    # pylab.plot([p1[0], p2[0], p3[0]], [p1[1], p2[1], p3[1]])

    x = []
    y = []

    if phi > psi:
        if fabs(phi - psi) > pi:
            psi += 2*pi
            r = arange(psi, phi, -0.001)
        else:
            r = arange(psi, phi, 0.001)
    else:
        if fabs(psi - phi) > pi:
            phi += 2*pi
            r = arange(psi, phi, 0.001)
        else:
            r = arange(psi, phi, -0.001)

    for a in r:
        x.append(c[0] + R*cos(a))
        y.append(c[1] + R*sin(a))

    # pylab.plot(x, y)

    return x, y


def LUR_R(V, B):
    return V**2/(g*tan(B))

def LUR(p1, p2, p3, R):
    angle1 = angle_line(p1, p2)
    angle2 = angle_line(p2, p3)

    turn_agnle = angle1 - angle2
    tas = sign(turn_agnle)
    if fabs(turn_agnle) > pi:
        turn_agnle -= tas*2*pi
    tas = sign(turn_agnle)

    lur = R * tan(fabs(turn_agnle)/2.0)


    n21 = normale(p2, p1)


    b = [
        p2[0] + lur*n21[0],
        p2[1] + lur*n21[1],
    ]

    c = [
        b[0] - tas*R*n21[1],
        b[1] + tas*R*n21[0],
    ]

    nb2 = normale(b, p2)
    nbc = normale(b, c)

    n23 = normale(p3, p2)
    
    d = moddot(n23, normale(c, p2))/mod(n23)

    e = [
        c[0] + tas*R*n23[1],
        c[1] - tas*R*n23[0],
    ]

    return b, e, c

def mod(v):
    return sqrt(v[0]**2 + v[1]**2)

def moddot(v, w):
    return fabs(v[0]*w[1] - v[1]*w[0])

def angle_line(p1, p2):
    return atan2(p2[1] - p1[1], p2[0] - p1[0])

def normale(p1, p2):
    n = [
        p2[0] - p1[0],
        p2[1] - p1[1]
    ]

    l = sqrt(n[0]**2 + n[1]**2)

    return [n[0] / l, n[1] / l]

if __name__ == "__main__":
    p = [
        [0, 0],
        [-100, -10],
        [-90, -100],
 
    ]
    V = 10

    pylab.figure()
    trajectory(p, V)
    pylab.show()