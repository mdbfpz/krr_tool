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

class TrajectoryPrediction:


    def run_simulation(self):#ovo mozda bolje da je static
        # Učitavanje podataka
        script_dir = os.path.dirname(os.path.abspath(__file__))

        ACsynonyms, FFP = load_mat_files() 

        # Simulation parameters
        Wind = [0, 0, 0]
        # SM = [10, 12.5, 15] # add safety margins

        brojevi_aviona=[1,2, 5, 10, 15, 20 , 25]
        vremena_sim=[1/60, 2/60, 5/60, 10/60, 15/60, 20/60, 25/60, 30/60] # hours
        
        for t_sim in vremena_sim:
            print(f"{int(t_sim * 60):6d}", end=" ")
        print("   minute",end=" ")
        for broj_aviona in brojevi_aviona:
            print()
            print(f"{broj_aviona:2d}", end=" ")
            for t_sim in vremena_sim:

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

                start = time.time()
                #for flight_name, flight_data in FFP.items():
                for i, (flight_name, flight_data) in enumerate(FFP.items()):
                    #print(flight_name)
                    if i == broj_aviona:  # Prekini nakon zeljenog broja aviona
                        break
                    #print(flight_name)
                    #print(json.dumps(flight_data,indent=4))
                    # Inicijalizacija načina rada zrakoplova
                    singleACarchive = np.zeros((SimulationTime, 21))
                    #print(f"velicina acarchive = {len(singleACarchive)}")
                    singleFlightACmode = {}
                    singleFlightACmode['ACC'] = 'C'  # Acceleration mode: (A)cceleration, (C)onstant, (D)ecceleration.
                    

                    #print(f"MOD PRIJE POSTAVLJANJA JE = {flight_data['mode']}")
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

                    def wrapper():
                        global singleACarchive, singleFlightACstate, singleFlightACcontrol, singleWPTi, singleFlightACmode
                        singleACarchive, singleFlightACstate, singleFlightACcontrol, singleWPTi, singleFlightACmode = trajectorygen_Weather_aware(
                            singleFlightACstate, singleFlightACcontrol, Wind, singleFlightACmode, dT, SimulationTime, WPTi, flight_data, 
                            FlightPath, opsdata, apfdata, GP, const, singleACarchive, 
                            endtime, desired_time
                        )
                    
                    
                    sim_start_time = time.time()
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

                    end = time.time()
                    
                    #print(str(i) + " " + str(broj_aviona))
                    if i == broj_aviona-1:  # Prekini nakon prvog (indeks 0)
                        print(f"{end - start:5.2f}s", end=" ")
                        pass
                        #break
                    """
                    print("\nOPSDATA:")
                    print(json.dumps(opsdata, indent=4, default=str))
                    print("=== REZULTATI SIMULACIJE ===")
                    
                    print(f"Flight: {flight_name}")
                    print("\nACmode:")
                    print(json.dumps(ACmode[flight_name], indent=4, default=str))
                    print("\nACcontrol:")
                    print(json.dumps(ACcontrol[flight_name], indent=4, default=str))
                    print("\nACstate:")
                    print(json.dumps(ACstate[flight_name], indent=4, default=str))
                    print("\nAtmHp:")
                    print(json.dumps(AtmHp[flight_name], indent=4, default=str))
                    print("\nACarchive (shape):")
                    print(f"  Shape: {ACarchive[flight_name].shape}")
                    print(f"  First 5 rows:\n{ACarchive[flight_name][:5]}")
                    print("\nWPTiDict:")
                    print(json.dumps(WPTiDict[flight_name], indent=4, default=str))
                    for key in ACarchive.keys():
                        print(key)
                    """
                    
"""
a = TrajectoryPrediction()
a.run_simulation()
"""
