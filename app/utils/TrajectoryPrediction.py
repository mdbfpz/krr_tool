# Load required data
import numpy as np
import scipy.io as sio
from math import radians, sin, cos, atan2, sqrt
import re
import os
import pprint
import json
from .trajPred_utils.trajectorygen import *
from .trajPred_utils.trajectorygen_utils import *
from .trajPred_utils.atmosphere import *
"""
from trajGen.trajPred_utils.trajectorygen import *
from trajGen.trajPred_utils.trajectorygen_utils import *
from trajGen.trajPred_utils.atmosphere import *
"""
import time
import cProfile


"""
mogu imat npr funkcije u ovoj klasi koje su slicne sa
load_mat_files al da čita s rdfoxa
preko tih datarepo dictionarya neceg

i vidit kaj onda s BADU podacima i kolko da ih ucitam i tak
kolko su informacije sinkronizirane po tipovima aviona i tak

mislim da sve podatke koje su prije bile iz .mat
mogu dobit iz rdfoxa za badu valjda opet bude badu nekakav samo ga
onda treba uskladit s ovim podacima
"""



"""
flight_data['mode']
        flight_data['type'] DA
        flight_data['mode'] spusta dize se stoji
        
        4d point:
        flight_data['xpos'],  # x koordinata (lonDeg) DA
        flight_data['ypos'],  # y koordinata (latDeg) DA
        flight_data['zpos'],  # visina (z koordinata) jos to fali
        flight_data['GS'],    # brzina (Ground Speed) DA
        flight_data['track'], # smjer (track angle) MOZDA
        ]FlightPath = flight_data['waypoints'] DA
"""
class TrajectoryPrediction:

    #constructor
    def __init__(self, lat, lon, VxMs, VyMs, icaoAircraftTypeDesignator, cleared_waypoints , current_fl  ):
        self.flight_data={}
        # Internally, x/y are expected later as lon/lat mapped to xpos/ypos
        # so keep explicit keys matching run_simulation's expectations.
        self.flight_data["xpos"]=lon  # lonDeg
        self.flight_data["ypos"]=lat  # latDeg
        self.flight_data["zpos"]=current_fl # visina (keep as provided; conversion handled elsewhere if needed)
        self.flight_data["GS"]=sqrt(VxMs**2 + VyMs**2)*1.94384 #pretvorba u cvorove
        self.flight_data["track"]=atan2(VxMs, VyMs)*180/np.pi #pretvorba u stupnjeve
        self.flight_data["type"]=icaoAircraftTypeDesignator
        self.flight_data["mode"]=2 #level
        self.flight_data["waypoints"]=self._parse_cleared_waypoints(cleared_waypoints)

    def _parse_cleared_waypoints(self, cleared_waypoints):
        """
        Normalize "cleared_waypoints" structure into internal FlightPath format.
        Output list items with keys: x (lon), y (lat), z (alt/FL, 0 if missing), name, flightoves=0, hist=1.

        Accepted inputs:
        - dict with structure like {"agreed": {"element": [ ... ]}}
        - direct list of elements with similar shape to items in agreed.element
        - JSON string path to file (optional convenience)
        """
        data = cleared_waypoints
        # If a string path is passed, try to load JSON
        if isinstance(data, str):
            try:
                if os.path.exists(data):
                    with open(data, "r", encoding="utf-8") as f:
                        data = json.load(f)
            except Exception:
                # Fall back to treating it as already-parsed
                pass

        # Drill down to elements array
        elements = []
        if isinstance(data, dict):
            if "agreed" in data and isinstance(data["agreed"], dict) and "element" in data["agreed"]:
                elements = data["agreed"]["element"] or []
            elif "element" in data:
                elements = data["element"] or []
        elif isinstance(data, list):
            elements = data

        flight_path = []
        for el in elements:
            try:
                # name
                name = None
                try:
                    name = el.get("elementStartPoint", {}).get("designatedPoint", {}).get("designator")
                except Exception:
                    name = None

                # position: pos string is typically "lat lon" (EPSG:4326)
                lat = None
                lon = None
                pos = (
                    el.get("point4D", {})
                      .get("position", {})
                      .get("pos")
                )
                if isinstance(pos, str):
                    parts = re.split(r"\s+", pos.strip())
                    if len(parts) >= 2:
                        # Convert to floats; assume order lat lon
                        try:
                            lat = float(parts[0])
                            lon = float(parts[1])
                        except ValueError:
                            lat = None
                            lon = None

                # level: flight level may be present; otherwise 0
                z = 0
                try:
                    lvl = (
                        el.get("point4D", {})
                          .get("level", {})
                          .get("flightLevel", {})
                          .get("value")
                    )
                    if lvl is not None:
                        # Keep as provided (e.g., FL). If conversion needed, handle downstream.
                        z = float(lvl)
                except Exception:
                    z = 0

                if lat is None or lon is None:
                    # Skip malformed entries
                    continue

                # Internal convention: x is lon, y is lat
                flight_path.append({
                    "x": lon,
                    "y": lat,
                    "z": z*100*0.3048,  # Convert FL to meters (1 FL = 100 feet; 1 foot = 0.3048 m)
                    "name": name,
                    "flightoves": 0,  # per requirement
                    "hist": 1         # per requirement
                })
            except Exception:
                # Ignore any malformed element
                continue

        return flight_path

    def run_simulation(self):#ovo mozda bolje da je static
        # Učitavanje podataka
        ACsynonyms = load_ACsynonyms() 

        FFP = load_FFS_data()

        # Simulation parameters
        Wind = [0, 0, 0]
        # SM = [10, 12.5, 15] # add safety margins

        t_sim=5/60#=[1/60, 2/60, 5/60, 10/60, 15/60, 20/60, 25/60, 30/60] # hours
        

        SimulationTime = int(t_sim * 3600 ) #int(1.5 * 3600 )
        desired_time = 7.75 * 3600  # start of simulation 7:45
        endtime = desired_time + SimulationTime  # end of simulation 9:15

        ACarchive = np.zeros((SimulationTime, 21))

        # all initial variables of aircraft modes must be defined on first
        # iteration and are set in ACmode struct.

        # Za svaki let u FFP
        ACmode = {}
        ACstate = {}
        ACcontrol = {}
        AtmHp = {}
        ACarchive = {}
        WPTiDict = {}


        flight_name, flight_data = [], []

        # Inicijalizacija načina rada zrakoplova
        singleACarchive = np.zeros((SimulationTime, 21))
        #print(f"velicina acarchive = {len(singleACarchive)}")
        singleFlightACmode = {}
        singleFlightACmode['ACC'] = 'C'  # Acceleration mode: (A)cceleration, (C)onstant, (D)ecceleration.
        

        # Postavljanje načina penjanja/spuštanja
        if flight_data['mode'] == 0:
            singleFlightACmode['CL'] = 'C'  # Climb mode: (C)limb
        elif flight_data['mode'] == 1:
            singleFlightACmode['CL'] = 'D'  # Descent mode: (D)escent
        elif flight_data['mode'] == 2:
            singleFlightACmode['CL'] = 'L'  # Level mode: (L)evel
        #print(f"MOD POSLIJE POSTAVLJANJA JE: {singleFlightACmode['CL']}")
        # Postavljanje ostalih parametara načina rada
        singleFlightACmode['ConfigMode'] = 'CL'  # Aircraft configuration mode: Clean
        singleFlightACmode['SpeedMode'] = 'C'   # Speed hold mode: CAS
        singleFlightACmode['tropo'] = False     # Is aircraft above tropopause?
        singleFlightACmode['incloud'] = 0       # It can be 0 (out) or 1 (in)
        
        # Početne varijable kretanja zrakoplova
        WPTi = 1  # Waypoint index. It marks the next waypoint in the waypoint list.
        singleFlightACmode['StillFlying'] = True
            
        # Pronalaženje BADA identifikatora za tip zrakoplova
        flight_type = flight_data['type']
        AC = None
        for j in range(len(ACsynonyms)):
            if ACsynonyms[j, 0] == flight_type:
                AC = ACsynonyms[j, 1]
                break
        
        # TODO: ne trebamo koristit sve A/C iz BADA-e, dobivat ćemo ih iz Oriona
        
        # Učitavanje BADA podataka iz OPF i APF datoteka
        # TODO: ovo maknuti iz petlje, loadati samo jednom i spremiti u varijablu dict i onda accessirati
        opsdata = read_OPF_file(f"{AC[0]}.OPF")  # Funkcija koju treba implementirati
        #print(opsdata)
        apfdata = read_APF_file(f"{AC[0]}.APF")  # Funkcija koju treba implementirati
        
        # Kreiranje početnog stanja zrakoplova
        #print(flight_data.keys())
        singleFlightACstate = [
            flight_data['xpos'],  # x koordinata (lonDeg)
            flight_data['ypos'],  # y koordinata (latDeg)
            flight_data['zpos'],  # visina (z koordinata)
            flight_data['GS'],    # brzina (Ground Speed)
            flight_data['track'], # smjer (track angle)
            opsdata['mref'] * 1000  # masa zrakoplova u kg
        ]

        # Provjera i prilagodba za male visine (vjerojatno polijetanje)
        if singleFlightACstate[2] < 150:
            singleFlightACmode['ACC'] = 'A'        # Ubrzavanje (Acceleration)
            singleFlightACmode['CL'] = 'C'         # Penjanje (Climb)
            singleFlightACmode['ConfigMode'] = 'T' # Konfiguracija za polijetanje (Take-off)
            singleFlightACstate[3] = 70            # Postavlja brzinu na 70 (vjerojatno čvorova)
            singleFlightACstate[2] = 1             # Postavlja visinu na 1 (vjerojatno metar ili stopu)
        
        ACmode[flight_name] = singleFlightACmode
        # Prilagodba mase za zrakoplove koji lete blizu svog plafona
        if singleFlightACstate[2] > opsdata['maxalt'] * 0.3048 * 0.94:
            singleFlightACstate[5] = singleFlightACstate[5] * 0.8
        
        # Dohvaćanje točaka rute i izračun maksimalne visine
        FlightPath = flight_data['waypoints']
        MaxAlt = max([wp['z'] for wp in FlightPath])
        
        # Prilagodba mase prema maksimalnoj visini i duljini leta
        if MaxAlt > opsdata['Hmax'] * 0.3048:
            singleFlightACstate[5] = 1000 * (opsdata['mref'] + opsdata['mmin']) / 2
        elif deg2nm(distance(FlightPath[0]['y'], FlightPath[0]['x'], FlightPath[-1]['y'], FlightPath[-1]['x'])) < 250:
            singleFlightACstate[5] = 1000 * ((opsdata['mref'] - opsdata['mmin']) * 0.2 + opsdata['mmin'])
        
        # Spremanje stanja za ovaj zrakoplov
        ACstate[flight_name] = singleFlightACstate
        """ovo je do 151 linije u matlab klasi main.m"""
        #print("\nACstate prije iceg:")
        #print(json.dumps(ACstate[flight_name], indent=4, default=str))
        singleFlightACcontrol = [300, 0, 0, 0]
        dT = 0

        singleFlightAtmHp = {}
        singleFlightAtmHp['T'], singleFlightAtmHp['p'], singleFlightAtmHp['rho'], singleFlightAtmHp['a'] = AtmosphereAtHp(singleFlightACstate[2], dT, const)
        CL = LiftCoefficient(singleFlightACstate[5], singleFlightAtmHp['rho'], opsdata['wingsurf'], singleFlightACstate[3], singleFlightACcontrol[1], const)
        CD = DragCoefficient(CL, opsdata['Cd0']['CR'], opsdata['Cd2']['CR'], opsdata['Cd0']['geardown'])
        
        singleFlightACcontrol[3] = CD
        ACcontrol[flight_name] = singleFlightACcontrol
        AtmHp[flight_name] = singleFlightAtmHp
        
        
        # simulating traffic
        if len(FlightPath) == 2:
            pass
        else:
            #cProfile.run('wrapper()')
            
            singleACarchive, singleFlightACstate, singleFlightACcontrol, singleWPTi, singleFlightACmode = trajectorygen_Weather_aware(
                singleFlightACstate, singleFlightACcontrol, Wind, singleFlightACmode, dT, SimulationTime, WPTi, flight_data, 
                FlightPath, opsdata, apfdata, GP, const, singleACarchive, 
                endtime, desired_time
            )
        #print("time spent inside simulation ", time.time()-sim_start_time)
        ACstate[flight_name] = singleFlightACstate   
        #print(f"singleFlightACcontrol[0] = {singleFlightACcontrol[0]}")     
        ACcontrol[flight_name] = singleFlightACcontrol
        #print(f"ACcontrol[flight_name] = {ACcontrol[flight_name]}")     
        AtmHp[flight_name] = singleFlightAtmHp
        ACarchive[flight_name] = singleACarchive
        WPTiDict[flight_name] = singleWPTi

            
                    
"""
a = TrajectoryPrediction()
a.run_simulation()
"""
