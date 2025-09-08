import numpy as np
from math import radians, sin, cos, atan2, sqrt
from .atmosphere import AirDensity,AtmosphereAtHp
import json
from .CDA import *





def deg2km(deg, radius=6371.0):
    """
    Pretvara kutnu udaljenost u stupnjevima u kilometre po velikoj kružnici na sferi.
    
    Parametri:
        deg : float ili array
            Kutna udaljenost u stupnjevima.
        radius : float, opcionalno
            Polumjer sfere (default: 6371.0 km, prosječni radijus Zemlje).
    
    Povrat:
        km : float ili array
            Udaljenost u kilometrima.
    """
    return np.asarray(deg) * (np.pi / 180.0) * radius

def ACCModeSet(TAS, desiredTAS, ACCmode, Vtol):
    """
    Određuje način akceleracije: (A)ccelerate, (C)onstant, (D)eccelerate.
    Sprječava trenutno prebacivanje s akceleracije na deakceleraciju i obrnuto.
    Ulazi:
        TAS        : trenutna stvarna brzina (m/s)
        desiredTAS : željena stvarna brzina (m/s)
        ACCmode    : trenutni način akceleracije ('A', 'C' ili 'D')
        Vtol       : tolerancija brzine (m/s)
    Izlaz:
        ACmode     : novi način akceleracije ('A', 'C' ili 'D')
    """
    if TAS is not None and desiredTAS is not None:
        if TAS < (desiredTAS - Vtol):  # Trenutna brzina manja od željene
            if ACCmode in ('C', 'A'):  # sprječava instantni prijelaz s deccel na accel
                return 'A'
            else:
                return 'C'
        elif TAS > (desiredTAS + Vtol):
            if ACCmode in ('C', 'D'):
                return 'D'
            else:
                return 'C'
        else:
            return 'C'

def tand(degrees):
    """
    Tangens od kuta u stupnjevima (ekvivalent MATLAB tand).
    """
    return np.tan(np.radians(degrees))


def trajectorygen_Weather_aware(
    ACstate, ACcontrol, Wind, ACmode, dT, SimulationTime, WPTi, FFP,
    waypoints, opsdata, apfdata, GP, const, ACarchive,
    endtime, desired_time, APlist=None
):
    """
    Generira putanju zrakoplova na temelju ulaznih parametara.
    """
    Loctime = desired_time

    # 1. Ako je maksimalna visina svih točaka rute < 160, vraća prazan arhiv
    if max([wp['z'] for wp in waypoints]) < 160:
        ACarchive = np.zeros((1, 21))
        return ACarchive, ACstate, ACcontrol, WPTi, ACmode

    WPtime = np.arange(1, SimulationTime + 2, 60)
    WPn = 1
    timestep = 1
    wp = 2
    ACalt = ACstate[2]
    ldg = 0
    cda = 0

    """
    # 2. Provjera je li zadnja točka rute aerodrom iz liste

    APlist_clean = []

    for item in APlist:
        if isinstance(item, (str, np.str_)):
            APlist_clean.append(str(item))
        elif isinstance(item, np.ndarray):
            if item.size > 0:
                APlist_clean.append(str(item.item()))  # extract the single value


    if len(FFP["waypoints"]) > 0:
        if FFP["waypoints"][-1]["name"] in APlist_clean:
            # Nije implementirano: CDA
            #print(f"accontrol[2] /////////////////////////7 ==== {ACcontrol[2]}")
            ACSm1 = ACstate.copy()
            ACcon1 = ACcontrol.copy()
            ACmode1 = ACmode.copy()
            waypoints = CDA(waypoints, ACSm1, WPTi, dT, const, ACmode1, GP, apfdata, opsdata, ACcon1, Wind)
            if waypoints is None:
                print("Clbo is zero, skipping simulation...")
                return None, None, None, None, None
            else:
                ldg = 1 """
                #print("kad se vratim iz cda u traj")
                #print(json.dumps(ACstate, indent=4, default=str))
                #print(f"accontrol[2] ==== {ACcontrol[2]}")
    # 3. Glavna simulacijska petlja
    #print("hačččččččč")
    inc = 0
    while ACmode['StillFlying']:
        inc += 1
        """ print("\nACstate prije iceg:")
        print(json.dumps(ACstate, indent=4, default=str)) """
        
        #print()
        #print(f"ACcontrol - potisak = {ACcontrol[0]}", flush=True)
        currenttime = Loctime + timestep - 1

        # 4. Atmosferski uvjeti na trenutnoj visini
        AtmHp_T, AtmHp_p, AtmHp_rho, AtmHp_a = AtmosphereAtHp(ACstate[2], dT, const)
        AtmHp = {'T': AtmHp_T, 'p': AtmHp_p, 'rho': AtmHp_rho, 'a': AtmHp_a}

        # 5. Provjera tropopauze
        ACmode['tropo'] = ACstate[2] > 11000

        # 6. Brzine
        CAScurrent = TAStoCAS(AtmHp['p'], AtmHp['rho'], ACstate[3])
        #print(ACstate[3])
        MACHcurrent = TAStoMach(AtmHp['T'], ACstate[3])
        #print(ACstate[3])
        TransAlt = CrossAlt(CAScurrent, MACHcurrent, dT, const)
        ACmode['SpeedMode'] = 'M' if ACstate[2] >= TransAlt else 'C'

        # 7. Faza leta i konfiguracija
        #print(len(ACstate))
        #ACmode['CL'] = CLModeSet(ACstate[2], waypoints[WPTi]['z'], ACmode['CL'], GP['Htol'])
        
        
        #print(f"TRAJECTORY GEN:::::ACstate[2]={ACstate[2]}, waypoints[WPTi]['z']={waypoints[WPTi]['z']}, ACmode['CL']={ACmode['CL']}, GP['Htol']={GP['Htol']}")

        #print(f"TRAJECTORY GEN::::: ACstate[2]={ACstate[2]}, CAScurrent={CAScurrent}, ACmode['CL']={ACmode['CL']}",flush=True)
        ACmode['CL'] = CLModeSet(ACstate[2], waypoints[WPTi]['z'], ACmode['CL'], GP['Htol'])
        

        ACmode['ConfigMode'] = ConfigModeSet(GP, opsdata, ACstate[2], CAScurrent, ACmode['CL'])
        
        #print(f"TRAJECTORY GEN:::::ACmode['ConfigMode'] = {ACmode['ConfigMode']}",flush=True)
        #print()
        # 8. Način ubrzanja i željena brzina
        desiredTAS = DesiredTASf(GP, apfdata, opsdata, ACmode['CL'], ACstate[5], AtmHp, ACmode['ConfigMode'], ACstate[2], TransAlt, const)
        if timestep == 1  and ACstate[3] > desiredTAS:
            #print(f"##########################################usporedba acstate[3] i desiredTAS: {ACstate[3]} i {desiredTAS}")
            ACstate[3] = desiredTAS
            #print(f"TG: ACstate[3]==={ACstate[3]}")
        ACmode['ACC'] = ACCModeSet(ACstate[3], desiredTAS, ACmode['ACC'], GP['Vtol'])

        # 9. Energy share factor
        esf = energysf(MACHcurrent, AtmHp['T'], dT, ACmode['ACC'], ACmode['CL'], ACmode['SpeedMode'], ACmode['tropo'], const)

        # 10. Kontrolni ulazi
        CL = LiftCoefficient(ACstate[5], AtmHp['rho'], opsdata['wingsurf'], ACstate[3], ACcontrol[1], const)
        ACcontrol[3] = ACDragCoef(opsdata, ACmode['ConfigMode'], CL)
        #print(f"prije 1. poziva ACStateUpdate acstate2=: {ACstate[2]}")
        ACcontrol[0] = ThrustSet(GP, opsdata, ACstate, ACcontrol, AtmHp, ACmode, desiredTAS, const)
        
        try:
            #print(f"ACcontrol prije poziva pitchset = {ACcontrol[2]}")
            ACcontrol[2] = PitchSet(opsdata, GP, ACstate, ACcontrol, waypoints[WPTi]['z'], esf, AtmHp, const)
            #print(f"ACcontrol pposlije poziva pitchset = {ACcontrol[2]}")
        except IndexError as e:
            #print(f"IndexError: {e}")
            #print(f"len(ACstate) = {len(ACcontrol)} (traženi indeks: 2)")
            #print(f"len(waypoints) = {len(waypoints)}, WPTi = {WPTi}")
            ACcontrol[2] = None  # ili neka zadana vrijednost

        
        ACcontrol[1] = BankSet(GP, ACstate, ACcontrol, waypoints, WPTi, ACmode['ConfigMode'], Wind, const)

        # 11. Potrošnja goriva
        
        FF = FuelConsumption(
            ACstate[3], ACcontrol[0], ACstate[2],
            opsdata['Cf1'], opsdata['Cf2'], opsdata['Cf3'], opsdata['Cf4'],
            opsdata['Cfcr'], opsdata['engtype'], ACmode['CL']
        )
        
        # 12. Arhiviranje stanja
        if desired_time <= currenttime <= endtime:
            ACarchive[timestep, 0] = ACstate[0]
            ACarchive[timestep, 1] = ACstate[1]
            ACarchive[timestep, 2] = ACstate[2]
            ACarchive[timestep, 3] = ACstate[3]
            ACarchive[timestep, 4] = ACstate[4]
            ACarchive[timestep, 5] = ACstate[5]
            ACarchive[timestep, 6] = ACcontrol[0]
            ACarchive[timestep, 7] = ACcontrol[1]
            ACarchive[timestep, 8] = ACcontrol[2]
            ACarchive[timestep, 9] = ACcontrol[3]
            ACarchive[timestep, 10] = desiredTAS
            ACarchive[timestep, 11] = esf
            ACarchive[timestep, 12] = CAScurrent
            ACarchive[timestep, 13] = MACHcurrent
            WindDir, WindSpd = cart2compass(-Wind[0], -Wind[1])
            track, GS = ACdrift(ACstate[4], ACstate[3], WindDir, WindSpd)
            ACarchive[timestep, 14] = track
            ACarchive[timestep, 15] = GS
            ACarchive[timestep, 16] = WPTi
            ACarchive[timestep, 17] = waypoints[WPTi]['x']
            ACarchive[timestep, 18] = waypoints[WPTi]['y']
            ACarchive[timestep, 19] = waypoints[WPTi]['z']
            ACarchive[timestep, 20] = currenttime

        # 13. Update ACstate
        #print(f"prije poziva ACStateUpdate: {ACstate[2]}")
        
        ACstate = ACStateUpdate(ACstate, ACcontrol, Wind, opsdata['wingsurf'], AtmHp['rho'], FF, CL, const, "tg")
        #print(f"nakon prvog poziva: ACstate[2]={ACstate[2]:.2f}, ACstate[3]={ACstate[3]:f},ACcontrol[2] = {ACcontrol[2]:f}, np.sin(ACcontrol[2])={np.sin(ACcontrol[2]):.4f}, Wind[2]={Wind[2]:.2f}, rezultat={ACstate[2] + ACstate[3] * np.sin(ACcontrol[2]) + Wind[2]:.2f}")
        
        #print(f"poslije poziva ACStateUpdate acstate2=: {ACstate[2]}")
        # 14. Provjera waypointa
        DistToNext = distance(ACstate[1], ACstate[0], waypoints[WPTi]['y'], waypoints[WPTi]['x'])
        DistToNext = deg2km(DistToNext) * 1000
        TurnRadius = ACstate[3] ** 2 / (const['g0'] * tand(GP['phi_nom_oth']))

        if waypoints[WPTi]['name'] == 'END':
            if distance(ACstate[1], ACstate[0], waypoints[WPTi]['y'], waypoints[WPTi]['x']) < 10:
                ACmode['StillFlying'] = False

        if waypoints[WPTi]['flyover'] == 1 and DistToNext < TurnRadius * 0.1:
            WPTi += 1
        elif waypoints[WPTi]['flyover'] == 0 and DistToNext < TurnRadius:
            WPTi += 1
        elif waypoints[WPTi]['name'] == 'addedWPT' and DistToNext < TurnRadius * 0.6:
            WPTi += 1

        # 15. Kraj simulacije?
        if WPTi >= len(waypoints):
            ACmode['StillFlying'] = False
        timestep += 1
        if timestep >= SimulationTime - 1:
            ACmode['StillFlying'] = False
        if currenttime >= endtime - 1:
            ACmode['StillFlying'] = False
        if ACstate[2] < -100:
            ACmode['StillFlying'] = False
        """ if inc == 5:
            break """
    #print(ACarchive, ACstate, ACcontrol, WPTi, ACmode)
    #print(f"na kraju je accontrol[0] = {ACcontrol[0]}")
    return ACarchive, ACstate, ACcontrol, WPTi, ACmode
