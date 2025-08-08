""" import scipy.io as sio
import numpy as np

# Učitavanje s opcijama za jednostavniji pristup
data = sio.loadmat('flight_pos.mat', struct_as_record=False, squeeze_me=True)
print(data)
flight_pos = data['flight_pos']

# Stvaranje rječnika za pohranu podataka o letovima
flights_data = {}

# Iteracija kroz sve letove i spremanje podataka
for flight in flight_pos:
    # Dohvaćanje imena leta
    flight_name = flight.name
    
    # Inicijalizacija liste za točke rute
    route_points = []
    
    # Dohvaćanje i spremanje svih točaka za trenutni let
    waypoints = flight.waypoints
    for wp in waypoints:
        point = {
            'x': wp.x,
            'y': wp.y,
            'z': wp.z,
            'flyover':wp.flyover,
            'hist':wp.hist
        }
        route_points.append(point)
    
    # Spremanje podataka o letu u glavni rječnik
    flights_data[flight_name] = {
        'waypoints': route_points,
        'type': flight.type if hasattr(flight, 'type') else None,
        'mode': flight.mode if hasattr(flight, 'mode') else None,
        'GS': flight.GS if hasattr(flight, 'GS') else None,
        'track': flight.track if hasattr(flight, 'track') else None
    }
#print(flights_data)
 # Primjer ispisa podataka za provjeru
for flight_name, flight_info in flights_data.items():
    print(f"Let: {flight_name}")
    print(f"Tip zrakoplova: {flight_info['type']}")
    print(f"Broj točaka rute: {len(flight_info['waypoints'])}")
    print("Prve tri točke rute:")
    for i, point in enumerate(flight_info['waypoints'][:3]):
        print(f"  Točka {i+1}: x={point['x']}, y={point['y']}, z={point['z']}, flyover={point['flyover']}")
    print("-" * 50) 
"""
# Direktan izračun
rezultat = (0.09963180809217054 * 122.6 * 1.224882290222165 * 70**2) / (2 * 51500.0)
#print(f"Rezultat: {rezultat:.6f}")



ACcontrol[0]=142300.96571539805, ACcontrol[2]=0.02177142857142857, ACcontrol[3]=0.09963106389408337
ACcontrol(1)=142300.96571539805, ACcontrol(3)=0.02177142857142857, ACcontrol(4)=0.09963104864322017


