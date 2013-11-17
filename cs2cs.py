#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
Computes UTM projection from WGS84 coordinates.
http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html
"""

from math import sin, cos, radians, sqrt, tan, floor

# Ellipsoid model constants (actual values for WGS84)
a = 6378137.0
b = 6356752.314
k0 = 0.9996  # UTMScaleFactor


def ArcLengthOfMeridian(phi):
    """
    Computes the ellipsoidal distance from the equator to a point at a
    given latitude.
   
    Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
    GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.
    
    Inputs:
         phi - Latitude of the point, in radians.
    Returns:
         The ellipsoidal distance of the point from the equator, in meters.
    """

    n = (a - b) / (a + b)

    alpha = (a + b) / 2.0 * (1.0 + n**2 / 4.0 + n**4 / 64.0)
    
    beta = -3.0 / 2.0 * n + 9.0 / 16.0 * n**3 - 3.0 / 32.0 * n**5
    
    gamma = 15.0 / 16.0 * n**2 - 15.0 / 32.0 * n**4
    
    delta = -35.0 / 48.0 * n**3 + 105.0 / 256.0 * n**5
    
    epsilon = 315.0 / 512.0 * n**4

    result = alpha * (phi + beta * sin(2.0 * phi) +
                        + gamma * sin(4.0 * phi) +
                        + delta * sin(6.0 * phi) +
                        + epsilon * sin(8.0 * phi))

    return result


def UTMCentralMeridian(zone):
    """
    Determines the central meridian for the given UTM zone.
    
    Inputs:
         zone - An integer value designating the UTM zone, range [1,60].
    
    Returns:
      The central meridian for the given UTM zone, in radians, or zero
      if the UTM zone parameter is outside the range [1,60].
      Range of the central meridian is the radian equivalent of [-177,+177].
    """
    cmeridian = radians(-183.0 + zone * 6.0)

    return cmeridian


def MapLatLonToXY(phi, lambda_, lambda0):
    """
    Converts a latitude/longitude pair to x and y coordinates in the
    Transverse Mercator projection.  Note that Transverse Mercator is not
    the same as UTM; a scale factor is required to convert between them.
    
    Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
    GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.
    
    Inputs:
       phi - Latitude of the point, in radians.
       lambda - Longitude of the point, in radians.
       lambda0 - Longitude of the central meridian to be used, in radians.

    Returns:    
       (x, y) - A 2-element tuple containing the x and y coordinates
                of the computed point.
    """
    ep2 = (a**2 - b**2) / b**2

    nu2 = ep2 * cos(phi)**2

    N = a**2 / (b * sqrt(1 + nu2))

    t = tan(phi)
    t2 = t * t

    l = lambda_ - lambda0

    l3 = 1.0 - t2 + nu2

    l4 = 5.0 - t2 + 9.0 * nu2 + 4.0 * nu2**2

    l5 = 5.0 - 18.0 * t2 + t2**2 + 14.0 * nu2 - 58.0 * t2 * nu2

    l6 = 61.0 - 58.0 * t2 + t2**2 + 270.0 * nu2 - 330.0 * t2 * nu2

    l7 = 61.0 - 479.0 * t2 + 179.0 * t2**2 - t2**3

    l8 = 1385.0 - 3111.0 * t2 + 543.0 * t2**2 - t2**3

    # Calculates easting x
    x = (N * cos(phi) * l + 
            N / 6.0 * cos(phi)**3 * l3 * l**3 + 
            N / 120.0 * cos(phi)**5 * l5 * l**5 + 
            N / 5040.0 * cos(phi)**7 * l7 * l**7)

    # Calculates nortinh y
    y = (ArcLengthOfMeridian(phi) +
            t / 2.0 * N * cos(phi)**2 * l**2 +
            t / 24.0 * N * cos(phi)**4 * l4 * l**4 +
            t / 720.0 * N * cos(phi)**6 * l6 * l**6 +
            t / 40320.0 * N * cos(phi)**8 * l * l**8)

    return (x, y)


def LatLonToUTMXY(lat, lon):
    """
    Converts a latitude/longitude pair to x and y coordinates in the
    Universal Transverse Mercator projection.
    
    Inputs:
       lat - Latitude of the point, in radians.
       lon - Longitude of the point, in radians.
    
    Returns:
       xy - A 2-element tuple where the UTM x and y values will be stored.
    """

    zone = floor((lon + 180.0) / 6.0) + 1
    x, y = MapLatLonToXY(lat, lon, UTMCentralMeridian(zone))

    x = x * k0 + 500000.0
    y = y * k0
    if y < 0.0:
        y += 10000000.0

    return (x, y)


if __name__ == "__main__":
    lat = radians(55.7522200)
    lon = radians(37.6155600)

    print(LatLonToUTMXY(lat, lon))
