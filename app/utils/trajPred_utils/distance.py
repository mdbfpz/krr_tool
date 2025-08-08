from math import radians, sin, cos, atan2, sqrt

def great_circle_deg(lat1, lon1, lat2, lon2):
    """
    Izračunava kutnu udaljenost (u stupnjevima) između dviju točaka na sferi.
    """
    # Pretvori u radijane
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    # Haversine formula
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    # Kutna udaljenost u radijanima, pretvori u stupnjeve
    deg = c * (180 / 3.141592653589793)
    return deg

def deg2nm(deg):
    """
    Pretvara stupnjeve luka u nautičke milje.
    """
    return deg * 60.0


deg_distance = great_circle_deg(40.71,-74.01,48.86,2.35)
#print(deg_distance)
nm_distance = deg2nm(deg_distance)
#print(nm_distance)
 
import numpy as np

def azimuth(lat1, lon1, lat2, lon2):
    """
    Vraća azimut (bearing) od točke (lat1, lon1) do (lat2, lon2), u radijanima.
    Ulazi i izlazi su u radijanima!
    Ako želiš stupnjeve, koristi np.degrees(az).
    """
    # Pretvori u radijane
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)

    dlon = lon2 - lon1

    x = np.sin(dlon) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)

    az = np.arctan2(x, y)  # rezultat u radijanima, od -pi do +pi
    az = (az + 2 * np.pi) % (2 * np.pi)  # normaliziraj na [0, 2pi)
    return az  # u radijanima

""" print(azimuth(10,10,10,40))
print(np.degrees(azimuth(10,10,10,40)))
 """

def reckon(lat1, lon1, arclen_deg, az_deg):
    """
    Računa novu točku (lat2, lon2) na sferi iz početne točke (lat1, lon1),
    duljine luka (arclen_deg, u stupnjevima) i azimuta (az_deg, u stupnjevima od sjevera CW).
    Povrat:
        lat2, lon2 : nova geografska širina i dužina (u stupnjevima)
    """
    # Pretvori ulaze u radijane
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    az = np.radians(az_deg)
    arclen = np.radians(arclen_deg)

    # Formula za reckoning na sferi (great-circle)
    lat2 = np.arcsin(np.sin(lat1) * np.cos(arclen) +
                     np.cos(lat1) * np.sin(arclen) * np.cos(az))
    lon2 = lon1 + np.arctan2(np.sin(az) * np.sin(arclen) * np.cos(lat1),
                             np.cos(arclen) - np.sin(lat1) * np.sin(lat2))

    # Pretvori natrag u stupnjeve
    lat2 = np.degrees(lat2)
    lon2 = np.degrees(lon2)
    return lat2, lon2

import math

def nm2deg(nm, sphere=None):
    """
    Convert distance from nautical miles to degrees
    
    Parameters:
    nm : float or array-like
        Distance in nautical miles
    sphere : float or str, optional
        Either radius in nautical miles or sphere name
        Default uses Earth's mean radius (3440.065 nm)
    
    Returns:
    float or array-like
        Distance in degrees
    """
    
    # Pozovi nm2rad funkciju
    if sphere is None:
        rad = nm2rad(nm)
    else:
        rad = nm2rad(nm, sphere)
    
    # Pretvori radijane u stupnjeve
    if hasattr(nm, '__iter__'):  # Ako je nm lista/array
        return [math.degrees(r) for r in rad]
    else:
        return math.degrees(rad)

def nm2rad(nm, sphere=None):
    """
    Convert distance from nautical miles to radians
    
    Parameters:
    nm : float or array-like
        Distance in nautical miles
    sphere : float or str, optional
        Either radius in nautical miles or sphere name
        Default uses Earth's mean radius (3440.065 nm)
    
    Returns:
    float or array-like
        Distance in radians
    """
    
    # Default Earth radius u nautičkim miljama
    earth_radius_nm = 3440.065
    
    if sphere is None:
        radius = earth_radius_nm
    elif isinstance(sphere, (int, float)):
        radius = sphere
    elif isinstance(sphere, str):
        # Sphere radijusi u nautičkim miljama
        sphere_radii = {
            'sun': 376006.11,
            'moon': 939.95,
            'mercury': 1516.0,
            'venus': 3760.4,
            'earth': 3440.065,
            'mars': 1830.0,
            'jupiter': 38792.0,
            'saturn': 32460.0,
            'uranus': 13724.0,
            'neptune': 13263.0,
            'pluto': 648.0
        }
        
        sphere_lower = sphere.lower()
        if sphere_lower in sphere_radii:
            radius = sphere_radii[sphere_lower]
        else:
            raise ValueError(f"Unknown sphere: {sphere}")
    else:
        raise ValueError("sphere must be a number or valid sphere name")
    
    # Pretvori nautičke milje u radijane
    if hasattr(nm, '__iter__'):  # Ako je nm lista/array
        return [n / radius for n in nm]
    else:
        return nm / radius



#print(nm2deg(500))
arclen = nm2deg(600)
az = 315
#print(reckon(51.5,0,arclen,az))
from math import pi
def tand(degrees):
    """
    Tangens od kuta u stupnjevima (ekvivalent MATLAB tand).
    """
    return np.tan(np.radians(degrees))

""" print(tand(90))   # 1.0
print(np.tan(np.pi/2))  """  # ogroman broj (npr. 1.633e+16, kao i u MATLAB-u)
