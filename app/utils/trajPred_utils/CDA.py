import numpy as np
from .atmosphere import AirDensity,AtmosphereAtHp
from math import radians, sin, cos, atan2, sqrt



def LiftCoefficient(mass, airdensity, wingarea, TAS, bankangle, const):
    """
    Izračunava koeficijent uzgona prema jednadžbi 3.6-1.
    
    mass       : masa zrakoplova [kg]
    airdensity : gustoća zraka [kg/m3]
    wingarea   : površina krila [m2]
    TAS        : stvarna brzina zraka [m/s]
    bankangle  : kut nagiba (bank angle) [radijani]
    const      : rječnik s konstantom g0
    """
    g0 = 9.806650
    CL = (2 * mass * g0) / (airdensity * TAS**2 * wingarea * np.cos(bankangle))
    #print(f'mass={mass:.6f}, g0={g0:}, airdensity={airdensity:.6f}, TAS={TAS:.6f}, wingarea={wingarea:.6f}, bankangle={bankangle:.6f}, CL={CL:.6f}')

    return CL

def distance(lat1, lon1, lat2, lon2):
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

GP = {
    'Almax': 2,
    'Anmax': 5,
    'phi_nom_to': 15,
    'phi_nom_oth': 30,
    'phi_nom_mil': 50,
    'phi_max_to': 25,
    'phi_max_hold': 35,
    'phi_max_oth': 45,
    'phi_max_mil': 70,
    'Cdes_exp': 1.6000,
    'Ctto': 1.2000,
    'Ctcr': 0.9500,
    'Hmax_to': 400,
    'Hmax_ic': 2000,
    'Hmax_ap': 8000,
    'Hmax_ld': 3000,
    'Cvmin_to': 1.2000,
    'Cvmin': 1.3000,
    'Vdcl1': 5,
    'Vdcl2': 10,
    'Vdcl3': 30,
    'Vdcl4': 60,
    'Vdcl5': 80,
    'Vdcl6': 20,
    'Vdcl7': 30,
    'Vdcl8': 35,
    'Vddes1': 5,
    'Vddes2': 10,
    'Vddes3': 20,
    'Vddes4': 50,
    'Vddes5': 5,
    'Vddes6': 10,
    'Vddes7': 20,
    'Vhold1': 230,
    'Vhold2': 240,
    'Vhold3': 265,
    'Vhold4': 0.8300,
    'Vbtrack': 35,
    'Vtaxi': 15,
    'Vapron': 10,
    'Vgate': 5,
    'Cred_tprop': 0.2500,
    'Cred_pis': 0,
    'Cred_jet': 0.1500,
    'Vtol': 1,
    'Htol': 15
}


def TAStoCAS(p, rho, TAS):
    """
    TAS to CAS conversion
    
    Inputs:
    p          :pressure in Pa.
    rho        :density in kg/m3
    TAS        :true airspeed m/s
    
    Output:
    cas        :calibrated airspeed m/s.
    
    REF: BADA 3.15      eq. 3.1-24
    """
    # Usklađeno s MATLAB verzijom - dodane zagrade oko cijele formule
    cas = (2*101325/(0.2857*1.225) * ((1 + (p/101325) * 
        ((1 + 0.2857*rho*TAS**2/(2*p))**3.5 - 1))**0.2857 - 1))**0.5
    return cas

def TAStoMach(T, TAS):
    """
    TAS to Mach conversion
    
    Inputs:
    T          :temp at alt in K
    TAS        :true airspeed m/s
    
    Output:
    Ma         :Mach number
    
    REF: BADA 3.15      eq. 3.1-26
    """
    # Usklađeno s MATLAB verzijom - koristi **0.5 umjesto np.sqrt
    Ma = TAS/((1.4*287.05287*T)**0.5)
    return Ma

def CrossAlt(CAS, Ma, dT, const):
    """
    calculates the altitude at which CAS expressed as TAS is equal to Mach
    number expressed as TAS. It is used to determine at which point in climb
    should the crossover from constant CAS climb to constant Mach climb occur.
    
    Inputs:
           CAS - calibrated airspeed, m/s
           Ma - Mach number
           dT - temperature deviation from ISA
           const - constants
    
    Output:
           Hcross - crossover altitude, m
    
    REF: BADA Manual 3.15. eq. 3.1-27 and 3.1-28
    """
    # Usklađeno s MATLAB verzijom - koristi **2 umjesto **2, dodane zagrade
    ptrans = const['P0']*((1+((const['kappa']-1)/2)*(CAS/const['a0'])**2)**(const['kappa']/(const['kappa']-1))-1)/ \
        ((1+((const['kappa']-1)*Ma**2)/2)**(const['kappa']/(const['kappa']-1))-1)

    Ttrop = const['T0']+dT+const['BetaT']*const['Hp_trop']

    ptrop = const['P0']*(((Ttrop-dT)/const['T0'])**(-const['g0']/(const['BetaT']*const['R'])))

    if ptrans>=ptrop:
       Hcross=(const['T0']/const['BetaT'])*((ptrans/const['P0'])**(-const['BetaT']*const['R']/const['g0'])-1)
    else:
       # Usklađeno s MATLAB verzijom - koristi np.log umjesto log
       Hcross=const['Hp_trop']-((const['R']*const['Tisa_trop'])/const['g0'])*np.log(ptrans/ptrop)

    return Hcross

def CLModeSet(h, h_desired, CLmode, Htol):
    """
    Determines the climb mode which can be (C)limb, (L)evel, and (D)escent.
    Prevents instantaneous switching from climb mode to descent and vice versa.
    input:
     h - Current geodetic altitude. [m]
     h_desired - Desired geodetic altitude. [m]
     CLmode - Current climb mode.
     Htol - tolerance from desired altitude in which mode will not be changed
     output:
    Climb mode which can be: (C)limb, (L)evel, and (D)escent.
    """
    if h < (h_desired - Htol):
        #print("OVO SE DOGODILOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO")# Current altitude is lower than desired altitude
        if CLmode == 'C' or CLmode == 'L':  # prevents instantly switching from climb to descent
            Cmode = 'C'
        else:
            Cmode = 'L'
    elif h > (h_desired + Htol):
        #print(f"MIJENJAM MODE U CLMODESET, TRENUTNI MOD JE {CLmode}",flush=True)
        if CLmode == 'D' or CLmode == 'L':
            Cmode = 'D'
        else:
            Cmode = 'L'
    else:
        Cmode = 'L'
    
    return Cmode


def ConfigModeSet(GP, opsdata, Hp, CAS, CLmode):
    """
    Određuje konfiguraciju zrakoplova: (T)ake-(O)ff, (I)nitial (C)limb, (CL)ean, (APP)roach, (L)an(D)in(G).
    """
    ConfigMode = 'ERR'
    VminCR = GP['Cvmin'] * opsdata['Vstall']['CR']
    VminAP = GP['Cvmin'] * opsdata['Vstall']['AP']

    if CLmode == 'C':
        if Hp < (GP['Hmax_to'] * 0.3048):
            ConfigMode = 'TO'
        elif (Hp >= (GP['Hmax_to'] * 0.3048)) and (Hp < (GP['Hmax_ic'] * 0.3048)):
            ConfigMode = 'IC'
        else:
            ConfigMode = 'CL'
    elif CLmode == 'L':
        ConfigMode = 'CL'
    else:
        if (Hp < (GP['Hmax_ld'] * 0.3048)) and (CAS < (VminAP + 10)):
            ConfigMode = 'LDG'
        elif (Hp < (GP['Hmax_ap'] * 0.3048)) and (Hp >= (GP['Hmax_ld'] * 0.3048)) and (CAS < (VminCR + 10)):
            ConfigMode = 'APP'
        elif (Hp <= (GP['Hmax_ld'] * 0.3048)) and (CAS < (VminCR + 10)) and (CAS >= (VminAP + 10)):
            ConfigMode = 'APP'
        elif (Hp >= (GP['Hmax_ap'] * 0.3048)):
            ConfigMode = 'CL'
        elif (Hp < (GP['Hmax_ap'] * 0.3048)) and (CAS >= (VminCR + 10)):
            ConfigMode = 'CL'
    return ConfigMode

def CAStoTAS(p, rho, CAS):
    """
    Pretvara Calibrated Air Speed (CAS) u True Air Speed (TAS).
    Ulazi:
        p   : tlak u Pa
        rho : gustoća zraka u kg/m3
        CAS : calibrated airspeed u m/s
    Izlaz:
        tas : true airspeed u m/s
    Referenca: BADA 3.15 eq. 3.1-23
    """
    tas = np.sqrt((2 * p / (0.2857 * rho)) * (
        (1 + (101325 / p) * ((1 + 0.2857 * 1.225 * CAS**2 / (2 * 101325))**3.5 - 1))**0.2857 - 1
    ))
    return tas


def SpeedforMass(Vref, mass, massref):
    """
    Izračun stvarne brzine za referentnu brzinu, referentnu masu i stvarnu masu (Eq. 3.4-1).
    Vref    : referentna brzina (bilo koja jedinica)
    mass    : stvarna masa zrakoplova (kg)
    massref : referentna masa zrakoplova (kg)
    Povrat: stvarna brzina (ista jedinica kao Vref)
    """
    return Vref * np.sqrt(mass / massref)

def MachtoTAS(T, Ma):
    """
    Pretvara Mach broj u True Air Speed (TAS).
    Ulazi:
        T  : temperatura na visini u K
        Ma : Mach broj
    Izlaz:
        TAS : true airspeed u m/s
    Referenca: BADA 3.15 eq. 3.1-26
    """
    TAS = Ma * np.sqrt(1.4 * 287.05287 * T)
    return TAS

def MinimumSpeed(Vstall, C_Vmin):
    """
    Izračunava minimalnu brzinu zrakoplova.
    Vstall  : brzina zastoja (stall speed) za fazu leta (bilo koja jedinica)
    C_Vmin  : koeficijent minimalne brzine (npr. 1.2 za TO, 1.3 inače)
    Povrat: minimalna brzina (ista jedinica kao Vstall)
    """
    return C_Vmin * Vstall

def JetLowSpeedBuffeting(BuffGradK, C_Lbo, weight, wingarea, pressure):
    """
    Računa granicu niske brzine buffetinga za jet zrakoplove.
    Vraća granicu kao Mach broj.
    """
    if BuffGradK == 0.0:
        BuffGradK = 0.1

    a1 = -C_Lbo / BuffGradK
    a3 = (weight / wingarea) / (0.583 * pressure * BuffGradK)

    """
    if(a1==0):
        print("o kurwa")
        exit(0)
    """
    Q = -(a1 ** 2) / 9

    R = (-27 * a3 - 2 * (a1 ** 3)) / 54

    # theta u radijanima
    theta = np.arccos(R / np.sqrt(-Q ** 3))

    # 120 i 240 stupnjeva u radijanima
    X1 = 2 * np.sqrt(-Q) * np.cos(theta / 3) - a1 / 3
    X2 = 2 * np.sqrt(-Q) * np.cos(theta / 3 + 2 * np.pi / 3) - a1 / 3
    X3 = 2 * np.sqrt(-Q) * np.cos(theta / 3 + 4 * np.pi / 3) - a1 / 3

    # Odabir najmanje pozitivne vrijednosti (kao u MATLAB-u)
    candidates = [X for X in [X1, X2, X3] if X > 0]
    if candidates:
        Mb = min(candidates)
    else:
        Mb = max(X1, X2, X3)  # ako su svi negativni, uzmi najveći (najmanje negativan)

    return Mb


def DesiredTASf(GP, apfdata, opsdata, CLmode, mass, AtmHp, ConfigMode, Hp, TransAlt, const):
    """
    Vraća željeni TAS prema BADA 3.15, Section 4.
    Indeksi prilagođeni MATLAB 1-based indeksaciji.
    """

    if CLmode == 'C':
        V_stall_ref = SpeedforMass(opsdata['Vstall']['TO'], mass, opsdata['mref']) / 0.5144

        if opsdata['engtype'].lower() in ['jet', 'j']:
            VCAS = [0]*8  # MATLAB VCAS(1) do VCAS(7) -> Python VCAS[0] do VCAS[6]
            VCAS[6] = apfdata['V_cl2']['AV']  # MATLAB VCAS(7)
            VCAS[5] = min(apfdata['V_cl1']['AV'], 250)  # MATLAB VCAS(6)
            VCAS[4] = min(GP['Cvmin'] * V_stall_ref + GP['Vdcl5'], VCAS[5])  # MATLAB VCAS(5)
            VCAS[3] = min(GP['Cvmin'] * V_stall_ref + GP['Vdcl4'], VCAS[4])  # MATLAB VCAS(4)
            VCAS[2] = min(GP['Cvmin'] * V_stall_ref + GP['Vdcl3'], VCAS[3])  # MATLAB VCAS(3)
            VCAS[1] = min(GP['Cvmin'] * V_stall_ref + GP['Vdcl2'], VCAS[2])  # MATLAB VCAS(2)
            VCAS[0] = min(GP['Cvmin'] * V_stall_ref + GP['Vdcl1'], VCAS[1])  # MATLAB VCAS(1)
            
            if Hp < 1500:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCAS[0] * 0.5144)  # MATLAB VCAS(1)
            elif Hp >= 1500 and Hp < 3000:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCAS[1] * 0.5144)  # MATLAB VCAS(2)
            elif Hp >= 3000 and Hp < 4000:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCAS[2] * 0.5144)  # MATLAB VCAS(3)
            elif Hp >= 4000 and Hp < 5000:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCAS[3] * 0.5144)  # MATLAB VCAS(4)
            elif Hp >= 5000 and Hp < 6000:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCAS[4] * 0.5144)  # MATLAB VCAS(5)
            elif Hp >= 6000 and Hp < 10000:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCAS[5] * 0.5144)  # MATLAB VCAS(6)
            elif Hp >= 10000 and Hp < TransAlt:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCAS[6] * 0.5144)  # MATLAB VCAS(7)
            elif Hp >= TransAlt:
                desiredTAS = MachtoTAS(AtmHp['T'], apfdata['M_cl']['AV'])

        elif opsdata['engtype'].lower() in ['turboprop', 't', 'piston', 'p']:
            VCASS = [0]*6  # MATLAB VCASS(1) do VCASS(5) -> Python VCASS[0] do VCASS[4]
            VCASS[4] = apfdata['V_cl2']['AV']  # MATLAB VCASS(5)
            VCASS[3] = min(apfdata['V_cl1']['AV'], 250)  # MATLAB VCASS(4)
            VCASS[2] = min(GP['Cvmin'] * V_stall_ref + GP['Vdcl8'], VCASS[3])  # MATLAB VCASS(3)
            VCASS[1] = min(GP['Cvmin'] * V_stall_ref + GP['Vdcl7'], VCASS[2])  # MATLAB VCASS(2)
            VCASS[0] = min(GP['Cvmin'] * V_stall_ref + GP['Vdcl6'], VCASS[1])  # MATLAB VCASS(1)
            
            if Hp < 500:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCASS[0] * 0.5144)  # MATLAB VCASS(1)
            elif Hp >= 500 and Hp < 1000:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCASS[1] * 0.5144)  # MATLAB VCASS(2)
            elif Hp >= 1000 and Hp < 1500:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCASS[2] * 0.5144)  # MATLAB VCASS(3)
            elif Hp >= 1500 and Hp < 10000:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCASS[3] * 0.5144)  # MATLAB VCASS(4)
            elif Hp >= 10000 and Hp < TransAlt:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCASS[4] * 0.5144)  # MATLAB VCASS(5)
            elif Hp >= TransAlt:
                desiredTAS = MachtoTAS(AtmHp['T'], apfdata['M_cl']['AV'])

    elif CLmode == 'L':
        if opsdata['engtype'].lower() in ['jet', 'j']:
            VCAS_ = [0]*5  # MATLAB VCAS_(1) do VCAS_(4) -> Python VCAS_[0] do VCAS_[3]
            VCAS_[3] = apfdata['V_cr2']['AV']  # MATLAB VCAS_(4)
            VCAS_[2] = min(apfdata['V_cr1']['AV'], 250)  # MATLAB VCAS_(3)
            VCAS_[1] = min(apfdata['V_cr1']['AV'], 220)  # MATLAB VCAS_(2)
            VCAS_[0] = min(apfdata['V_cr1']['AV'], 170)  # MATLAB VCAS_(1)
            
            if Hp < 3000:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCAS_[0] * 0.5144)  # MATLAB VCAS_(1)
            elif Hp >= 3000 and Hp < 6000:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCAS_[1] * 0.5144)  # MATLAB VCAS_(2)
            elif Hp >= 6000 and Hp < 14000:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCAS_[2] * 0.5144)  # MATLAB VCAS_(3)
            elif Hp >= 14000 and Hp < TransAlt:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCAS_[3] * 0.5144)  # MATLAB VCAS_(4)
            elif Hp >= TransAlt:
                desiredTAS = MachtoTAS(AtmHp['T'], apfdata['M_cr']['AV'])

        elif opsdata['engtype'].lower() in ['turboprop', 't', 'piston', 'p']:
            VCASa = [0]*5  # MATLAB VCASa(1) do VCASa(4) -> Python VCASa[0] do VCASa[3]
            VCASa[3] = apfdata['V_cr2']['AV']  # MATLAB VCASa(4)
            VCASa[2] = min(apfdata['V_cr1']['AV'], 250)  # MATLAB VCASa(3)
            VCASa[1] = min(apfdata['V_cr1']['AV'], 180)  # MATLAB VCASa(2)
            VCASa[0] = min(apfdata['V_cr1']['AV'], 150)  # MATLAB VCASa(1)
            
            if Hp < 3000:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCASa[0] * 0.5144)  # MATLAB VCASa(1)
            elif Hp >= 3000 and Hp < 6000:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCASa[1] * 0.5144)  # MATLAB VCASa(2)
            elif Hp >= 6000 and Hp < 10000:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCASa[2] * 0.5144)  # MATLAB VCASa(3)
            elif Hp >= 10000 and Hp < TransAlt:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCASa[3] * 0.5144)  # MATLAB VCASa(4)
            elif Hp >= TransAlt:
                desiredTAS = MachtoTAS(AtmHp['T'], apfdata['M_cr']['AV'])
        #print(f"DESIREDTAS = {desiredTAS}")
    elif CLmode == 'D':
        V_stall_ref = SpeedforMass(opsdata['Vstall']['LD'], mass, opsdata['mref']) / 0.5144
        
        if opsdata['engtype'].lower() in ['jet', 'j']:
            VCAS = [0]*8  # MATLAB VCAS(1) do VCAS(7) -> Python VCAS[0] do VCAS[6]
            VCAS[6] = apfdata['V_des2']['AV']  # MATLAB VCAS(7)
            VCAS[5] = min(apfdata['V_des1']['AV'], 250)  # MATLAB VCAS(6)
            VCAS[4] = min(apfdata['V_des1']['AV'], 220)  # MATLAB VCAS(5)
            VCAS[3] = min(GP['Cvmin'] * V_stall_ref + GP['Vddes4'], VCAS[4])  # MATLAB VCAS(4)
            VCAS[2] = min(GP['Cvmin'] * V_stall_ref + GP['Vddes3'], VCAS[3])  # MATLAB VCAS(3)
            VCAS[1] = min(GP['Cvmin'] * V_stall_ref + GP['Vddes2'], VCAS[2])  # MATLAB VCAS(2)
            VCAS[0] = min(GP['Cvmin'] * V_stall_ref + GP['Vddes1'], VCAS[1])  # MATLAB VCAS(1)
            
            if Hp < 1000:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCAS[0] * 0.5144)  # MATLAB VCAS(1)
            elif Hp >= 1000 and Hp < 1500:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCAS[1] * 0.5144)  # MATLAB VCAS(2)
            elif Hp >= 1500 and Hp < 2000:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCAS[2] * 0.5144)  # MATLAB VCAS(3)
            elif Hp >= 2000 and Hp < 3000:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCAS[3] * 0.5144)  # MATLAB VCAS(4)
            elif Hp >= 3000 and Hp < 6000:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCAS[4] * 0.5144)  # MATLAB VCAS(5)
            elif Hp >= 6000 and Hp < 10000:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCAS[5] * 0.5144)  # MATLAB VCAS(6)
            elif Hp >= 10000 and Hp < TransAlt:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCAS[6] * 0.5144)  # MATLAB VCAS(7)
            elif Hp >= TransAlt:
                desiredTAS = MachtoTAS(AtmHp['T'], apfdata['M_des']['AV'])

        elif opsdata['engtype'].lower() in ['turboprop', 't', 'piston', 'p']:
            VCASS = [0]*6  # MATLAB VCASS(1) do VCASS(5) -> Python VCASS[0] do VCASS[4]
            VCASS[4] = apfdata['V_des2']['AV']  # MATLAB VCASS(5)
            VCASS[3] = apfdata['V_des1']['AV']  # MATLAB VCASS(4)
            VCASS[2] = min(GP['Cvmin'] * V_stall_ref + GP['Vddes7'], VCASS[3])  # MATLAB VCASS(3)
            VCASS[1] = min(GP['Cvmin'] * V_stall_ref + GP['Vddes6'], VCASS[2])  # MATLAB VCASS(2)
            VCASS[0] = min(GP['Cvmin'] * V_stall_ref + GP['Vddes5'], VCASS[1])  # MATLAB VCASS(1)
            
            if Hp < 500:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCASS[0] * 0.5144)  # MATLAB VCASS(1)
            elif Hp >= 500 and Hp < 1000:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCASS[1] * 0.5144)  # MATLAB VCASS(2)
            elif Hp >= 1000 and Hp < 1500:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCASS[2] * 0.5144)  # MATLAB VCASS(3)
            elif Hp >= 1500 and Hp < 10000:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCASS[3] * 0.5144)  # MATLAB VCASS(4)
            elif Hp >= 10000 and Hp < TransAlt:
                desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], VCASS[4] * 0.5144)  # MATLAB VCASS(5)
            elif Hp >= TransAlt:
                desiredTAS = MachtoTAS(AtmHp['T'], apfdata['M_des']['AV'])

    # Provjera minimalne i maksimalne dopuštene brzine
    if ConfigMode == 'TO':
        V_stall_ref = SpeedforMass(opsdata['Vstall']['TO'], mass, opsdata['mref'])
        SpeedMin = MinimumSpeed(V_stall_ref, GP['Cvmin_to'])
    elif ConfigMode == 'IC':
        V_stall_ref = SpeedforMass(opsdata['Vstall']['IC'], mass, opsdata['mref'])
        SpeedMin = MinimumSpeed(V_stall_ref, GP['Cvmin'])
    elif ConfigMode == 'CL':
        V_stall_ref = SpeedforMass(opsdata['Vstall']['CR'], mass, opsdata['mref'])
        SpeedMin = MinimumSpeed(V_stall_ref, GP['Cvmin'])
    elif ConfigMode == 'APP':
        V_stall_ref = SpeedforMass(opsdata['Vstall']['AP'], mass, opsdata['mref'])
        SpeedMin = MinimumSpeed(V_stall_ref, GP['Cvmin'])
    elif ConfigMode == 'LDG':
        V_stall_ref = SpeedforMass(opsdata['Vstall']['LD'], mass, opsdata['mref'])
        SpeedMin = MinimumSpeed(V_stall_ref, GP['Cvmin'])

    # Low speed buffeting limit
    if Hp > (15000 * 0.3048) and opsdata['engtype'] == 'Jet':
        """
        if(opsdata['Clbo']==0):
            print("o kurwa")
            exit(0)
        """
        BuffetingLimitMach = JetLowSpeedBuffeting(
            opsdata['k'], opsdata['Clbo'], const['g0'] * mass,
            opsdata['wingsurf'], AtmHp['p']
        )
        if SpeedMin < MachtoTAS(AtmHp['T'], BuffetingLimitMach):
            SpeedMin = MachtoTAS(AtmHp['T'], BuffetingLimitMach)

    # Finalne provjere
    if desiredTAS < SpeedMin:
        desiredTAS = SpeedMin

    if desiredTAS > CAStoTAS(AtmHp['p'], AtmHp['rho'], opsdata['VMO'] * 0.5144):
        desiredTAS = CAStoTAS(AtmHp['p'], AtmHp['rho'], opsdata['VMO'] * 0.5144)

    if desiredTAS > MachtoTAS(AtmHp['T'], opsdata['MMO']):
        desiredTAS = MachtoTAS(AtmHp['T'], opsdata['MMO'])

    return desiredTAS


def ACCModeSet(TAS, desiredTAS, ACCmode, Vtol):
    """
    Određuje način akceleracije: (A)ccelerate, (C)onstant, (D)eccelerate.
    Sprječava trenutno prebacivanje s akceleracije na deakceleraciju i obrnuto.

    Parametri:
        TAS        : trenutna stvarna brzina (m/s)
        desiredTAS : željena stvarna brzina (m/s)
        ACCmode    : trenutni način akceleracije ('A', 'C' ili 'D')
        Vtol       : tolerancija brzine (m/s)
    Povrat:
        ACmode     : novi način akceleracije ('A', 'C' ili 'D')
    """
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

def energysf(Ma, T, dt, ACCmode, CLmode, SpeedMode, Tropopause, const):
    """
    Calculates Energy Share Factor as a function of Mach number and climb/acceleration mode.
    Inputs:
        Ma          - Mach number. [dimensionless]
        T           - Air temperature. [K]
        dt          - Temperature difference from ISA conditions. [K]
        ACCmode     - Acceleration mode: (A)ccelerating, (C)onstant, or (D)ecelerating.
        CLmode      - Climb mode: (C)limb, (L)evel, (D)escent
        SpeedMode   - Speed mode: (C)AS, (M)ach
        Tropopause  - Is aircraft above tropopause? (True/False)
        const       - constants object or dict with kappa, R, BetaT, g0
    Output:
        esf         - energy share factor
    REF: BADA User Manual 3.15 eq. 3.2-8 to 3.2-11
    """
    # Podrška za dict i objekt
    kappa = const['kappa'] if isinstance(const, dict) else const.kappa
    R = const['R'] if isinstance(const, dict) else const.R
    BetaT = const['BetaT'] if isinstance(const, dict) else const.BetaT
    g0 = const['g0'] if isinstance(const, dict) else const.g0

    a = (kappa * R * BetaT * Ma * Ma * (T - dt)) / (2 * g0 * T)
    b = 1 + (kappa - 1) * Ma * Ma / 2

    if CLmode == 'L':  # Level
        esf = 0
    elif CLmode == 'C':  # Climb
        if ACCmode == 'C':  # Constant speed
            if SpeedMode == 'C' and Tropopause:
                esf = 1 / (1 + b ** (-1 / (kappa - 1)) * (b ** (kappa / (kappa - 1)) - 1))
            elif SpeedMode == 'C' and not Tropopause:
                esf = 1 / (1 + a + b ** (-1 / (kappa - 1)) * (b ** (kappa / (kappa - 1)) - 1))
            elif SpeedMode == 'M' and Tropopause:
                esf = 1
            elif SpeedMode == 'M' and not Tropopause:
                esf = 1 / (1 + a)
            else:
                esf = None
        elif ACCmode == 'A':  # Acceleration
            esf = 0.5  # ili 0.3 prema komentaru
        else:  # Deceleration
            esf = 1.5  # ili 1.7 prema komentaru
    else:  # Descent
        if ACCmode == 'C':  # Constant speed
            if SpeedMode == 'C' and Tropopause:
                esf = 1 / (1 + b ** (-1 / (kappa - 1)) * (b ** (kappa / (kappa - 1)) - 1))
            elif SpeedMode == 'C' and not Tropopause:
                esf = 1 / (1 + a + b ** (-1 / (kappa - 1)) * (b ** (kappa / (kappa - 1)) - 1))
            elif SpeedMode == 'M' and Tropopause:
                esf = 1
            elif SpeedMode == 'M' and not Tropopause:
                esf = 1 / (1 + a)
            else:
                esf = None
        elif ACCmode == 'A':  # Acceleration
            esf = 1.5  # ili 1.7 prema komentaru
        else:  # Deceleration
            esf = 0.5  # ili 0.3 prema komentaru
    return esf

def DragCoefficient(liftcoefficient, C_D0, C_D2, C_D0dLDG):
    """
    Izračunava koeficijent otpora (CD).
    """
    #print(f'C_D0={C_D0:.6f}, C_D0dLDG={C_D0dLDG:.6f}, C_D2={C_D2:.6f}, liftcoefficient={liftcoefficient:.6f}, rezultat={C_D0 + C_D0dLDG + C_D2 * liftcoefficient ** 2:.6f}')

    return C_D0 + C_D0dLDG + C_D2 * liftcoefficient ** 2


def ACDragCoef(opsdata, ConfMode, CL):
    """
    Određuje koeficijent otpora (CD) na temelju konfiguracije i koeficijenta uzgona.
    """
    # Pretpostavlja se da su opsdata['Cd0'] i opsdata['Cd2'] rječnici s ključevima 'AP', 'LD', 'CR', itd.
    if ConfMode == 'APP' and opsdata['Cd0']['AP'] != 0:
        #print("prvi")
        DragCoef = DragCoefficient(CL, opsdata['Cd0']['AP'], opsdata['Cd2']['AP'], 0)
    elif ConfMode == 'LDG' and opsdata['Cd0']['LD'] != 0:
        #print("drugi")
        DragCoef = DragCoefficient(CL, opsdata['Cd0']['LD'], opsdata['Cd2']['LD'], opsdata['Cd0']['geardown'])
    else:
        #print("treci")
        DragCoef = DragCoefficient(CL, opsdata['Cd0']['CR'], opsdata['Cd2']['CR'], 0)
    return DragCoef

def ThrustMaxClimb(Hp, dT, C_Tc1, C_Tc2, C_Tc3, C_Tc4, C_Tc5, TAS, EngType):
    """
    Računa maksimalni potisak u penjanju (Tmaxcl) prema BADA.
    """
    Hp_ft = Hp / 0.3048
    TAS_kt = TAS * 3600 / 1852

    EngType = EngType.lower()
    if EngType in ['jet', 'j']:
        TmaxclISA = C_Tc1 * (1 - Hp_ft / C_Tc2 + C_Tc3 * Hp_ft ** 2)
    elif EngType in ['turboprop', 't']:
        TmaxclISA = (C_Tc1 / TAS_kt) * (1 - Hp_ft / C_Tc2) + C_Tc3
    elif EngType in ['piston', 'p']:
        TmaxclISA = C_Tc1 * (1 - Hp_ft / C_Tc2) + C_Tc3 / TAS_kt
    else:
        TmaxclISA = -99999999999999

    dTeff = dT - C_Tc4
    if C_Tc5 < 0:
        C_Tc5 = 0
    a = dTeff * C_Tc5
    if a < 0:
        a = 0
    if a > 0.4:
        a = 0.4

    Tmaxcl = TmaxclISA * (1 - a)
    return Tmaxcl

def ThrustInDescent(Hp, Hpdes, Tmaxcl, C_Tdeslow, C_Tdeshigh, C_Tdesapp, C_Tdesld, Configuration):
    """
    Računa potisak u spuštanju (Thrust in Descent).
    """
    #print()
    #print(Hp,Hpdes)
    #print("HOKI KE BUUUUUUUUUUUUUUUUUUU")
    """ if abs(Hp - Hpdes) > 5:
        print("HAAAAAAAAAAAAAAAAAAAAAAA")
        TiD = C_Tdeshigh * Tmaxcl """
    if Hp - Hpdes:
        #print("HAAAAAAAAAAAAAAAAAAAAAAA")
        TiD = C_Tdeshigh * Tmaxcl
    else:
        conf = Configuration.upper()
        if conf in ['CL', 'CLEAN']:
            #print("hahaha")
            #print(C_Tdeslow,Tmaxcl)
            TiD = C_Tdeslow * Tmaxcl
        elif conf in ['APP', 'APPROACH']:
            TiD = C_Tdesapp * Tmaxcl
        elif conf in ['LDG', 'LANDING']:
            TiD = C_Tdesld * Tmaxcl
        else:
            TiD = -999999999
    return TiD
#print(ThrustInDescent(29899.70,29831.00,67203.37,0.108470,0.136030,0.157490,0.395660,"CL"))
def ThrustMaxCruise(Tmaxcl, Ctcr):
    """
    Računa maksimalni mogući potisak u krstarenju.
    """
    return Tmaxcl * Ctcr

def MaximumAltitude(h_MO, hmax, Gt, Gw, dT, C_tc4, mass_max, mass_actual):
    """
    Izračunava maksimalnu visinu za zadanu masu i temperaturu.
    Svi parametri su u istim jedinicama kao u MATLAB-u.
    Povrat: maksimalna visina u stopama.
    """
    a = dT - C_tc4
    if a < 0:
        a = 0

    b = hmax + Gt * a + Gw * (mass_max - mass_actual)

    if hmax == 0:
        maxalt = h_MO  # prema BADA 3.15, poglavlje 3.5.
    else:
        if h_MO < b:
            maxalt = h_MO
        else:
            maxalt = b
    return maxalt

def ReducedClimbPowerCoeff(hmax, Hp, m_max, m_min, m_act, ACtype, GP):
    """
    Izračunava koeficijent smanjene snage penjanja prema BADA 3.15, Eq. 3.8-2.
    
    Parametri:
        hmax   : maksimalna visina (feet)
        Hp     : trenutna geopotencijalna visina (feet)
        m_max  : maksimalna masa (kg)
        m_min  : minimalna masa (kg)
        m_act  : trenutna masa (kg)
        ACtype : tip motora ('P', 'Piston', 'T', 'Turboprop', 'J', 'Jet')
        GP     : rječnik ili objekt s članovima Cred_pis, Cred_tprop, Cred_jet
    Povrat:
        Cpowered : koeficijent smanjene snage penjanja (float)
    """
    if Hp < hmax * 0.8:
        ac_type = ACtype.upper()
        if ac_type in ['P', 'PISTON']:
            Cred = GP['Cred_pis'] if isinstance(GP, dict) else GP.Cred_pis
        elif ac_type in ['T', 'TURBOPROP']:
            Cred = GP['Cred_tprop'] if isinstance(GP, dict) else GP.Cred_tprop
        elif ac_type in ['J', 'JET']:
            Cred = GP['Cred_jet'] if isinstance(GP, dict) else GP.Cred_jet
        else:
            Cred = -999999999
    else:
        Cred = 0

    Cpowered = 1 - Cred * (m_max - m_act) / (m_max - m_min)
    return Cpowered


def ThrustSet(GP, opsdata, ACstate, ACcontrol, AtmHp, ACmode, DesTAS, const):
    """
    Određuje potrebni potisak (thrust) u Newtonima.
    """
    Thrust = -999999
    ACcontrol[3] = round(ACcontrol[3], 6)
    TmaxCl = ThrustMaxClimb(
        ACstate[2], AtmHp['T'] - 288.15, opsdata['Ctc1'], opsdata['Ctc2'],
        opsdata['Ctc3'], opsdata['Ctc4'], opsdata['Ctc5'], ACstate[3], opsdata['engtype']
    )
    #print(f"ThrustInDescent args: ACstate[2]/0.3048={ACstate[2]/0.3048}, Hpdes={opsdata['Hpdes']}, TmaxCl={TmaxCl}, Ctdeslow={opsdata['Ctdeslow']}, Ctdeshi={opsdata['Ctdeshi']}, Ctdesapp={opsdata['Ctdesapp']}, Ctdesld={opsdata['Ctdesld']}, ConfigMode={ACmode['ConfigMode']}")

    Tdes = ThrustInDescent(
        ACstate[2] / 0.3048, opsdata['Hpdes'], TmaxCl, opsdata['Ctdeslow'],
        opsdata['Ctdeshi'], opsdata['Ctdesapp'], opsdata['Ctdesld'], ACmode['ConfigMode']
    )
    #print(f"TDES={Tdes}")
    TmaxCruise = ThrustMaxCruise(TmaxCl, GP['Ctcr'])

    Hmax = MaximumAltitude(
        opsdata['maxalt'], opsdata['Hmax'], opsdata['tempgrad'], opsdata['gw'],
        AtmHp['T'] - 288.15, opsdata['Ctc4'], opsdata['mmax'], ACstate[5]
    )
    CTred = ReducedClimbPowerCoeff(
        Hmax, ACstate[2], opsdata['mmax'], opsdata['mmin'], ACstate[5], opsdata['engtype'], GP
    )

    if ACmode['CL'] == 'L' and ACmode['ACC'] == 'C':
        #print(1)
        ThrForSpeed = (ACcontrol[3] * opsdata['wingsurf'] * AtmHp['rho'] * DesTAS ** 2) / 2 + const['g0'] * np.sin(ACcontrol[2]) / ACstate[5]
        Thrust = ThrForSpeed
        if ThrForSpeed > TmaxCruise:
            Thrust = TmaxCruise

    elif ACmode['CL'] == 'L' and ACmode['ACC'] == 'A':
        #print(2)
        ThrForMaxAccel = GP['Almax'] * 0.3048 * ACstate[5] + (ACcontrol[3] * opsdata['wingsurf'] * AtmHp['rho'] * ACstate[3] ** 2) / 2 + const['g0'] * np.sin(ACcontrol[2]) * ACstate[5]
        Thrust = ThrForMaxAccel
        if ThrForMaxAccel > TmaxCruise:
            Thrust = TmaxCruise

    elif ACmode['CL'] == 'L' and ACmode['ACC'] == 'D':
        #print(3)
        ThrForMaxDecel = -GP['Almax'] * 0.3048 * ACstate[5] + (ACcontrol[3] * opsdata['wingsurf'] * AtmHp['rho'] * ACstate[3] ** 2) / 2 + const['g0'] * np.sin(ACcontrol[2]) * ACstate[5]
        Thrust = ThrForMaxDecel
        if ThrForMaxDecel > TmaxCruise:
            Thrust = TmaxCruise

    elif ACmode['CL'] == 'C':
        #print(4)
        if ACstate[2] < (Hmax * 0.8 * 0.3048):
            Thrust = TmaxCl * CTred
        else:
            Thrust = TmaxCl
        if Thrust > TmaxCl:
            Thrust = TmaxCl

    elif ACmode['CL'] == 'D':
        #print(5)
        Thrust = Tdes
        if Thrust < 0:
            Thrust = 0
        #print(Thrust)
    #print(Thrust,flush=True)
    return Thrust


def PitchSet(opsdata, GP, ACstate, ACcontrol, DesiredAlt, ESF, AtmHp, const):
    """
    Determines the required pitch angle (in radians), usklađeno s MATLAB logikom.
    """
    CurrentAlt = ACstate[2]
    CurrentSpd = ACstate[3]
    CurrentPitch = ACcontrol[2]
    Thr = ACcontrol[0]
    CD = ACcontrol[3]
    Mass = ACstate[5]

    dPitchMax = GP['Anmax'] * 0.3048 / CurrentSpd

    AltToLevel = 0
    NSteps = CurrentPitch / dPitchMax
    NSteps = abs(round(NSteps))
    dAlt = abs(DesiredAlt - CurrentAlt)

    if NSteps == 0:
        AltToLevel = GP['Htol']
    else:
        i = 0
        while i <= NSteps:
            AltToLevel += CurrentSpd * np.sin(abs(abs(CurrentPitch) - i * dPitchMax))
            i += 1

    # Izračun DesiredPitch (ograniči asin na [-1,1])
    val = (Thr - (CD * opsdata['wingsurf'] * AtmHp['rho'] * CurrentSpd ** 2) / 2) * ESF / (const['g0'] * Mass)
    val = np.clip(val, -1, 1)
    DesiredPitch = np.arcsin(val)

    # === Ovdje je ključ: nemoj dozvoliti negativan pitch ako si ispod željene visine ===
    if CurrentAlt < (DesiredAlt - GP['Htol']):
        if dAlt <= AltToLevel:
            if CurrentPitch >= 0:
                NewPitch = max(0, CurrentPitch - dPitchMax)
            else:
                NewPitch = min(0, CurrentPitch + dPitchMax)
        else:
            if CurrentPitch < (DesiredPitch - dPitchMax):
                NewPitch = CurrentPitch + dPitchMax
            elif CurrentPitch > (DesiredPitch + dPitchMax):
                NewPitch = CurrentPitch - dPitchMax
            else:
                NewPitch = max(0, DesiredPitch)  # Ovdje forsiramo samo pozitivne vrijednosti za climb
    elif CurrentAlt > (DesiredAlt + GP['Htol']):
        if dAlt <= AltToLevel:
            if CurrentPitch <= 0:
                NewPitch = CurrentPitch + dPitchMax
            else:
                NewPitch = CurrentPitch - dPitchMax
        else:
            if CurrentPitch < (DesiredPitch - dPitchMax):
                NewPitch = CurrentPitch + dPitchMax
            elif CurrentPitch > (DesiredPitch + dPitchMax):
                NewPitch = CurrentPitch - dPitchMax
            else:
                NewPitch = min(0, DesiredPitch)  # Ovdje forsiramo samo negativne vrijednosti za descent
    else:
        if CurrentPitch < -dPitchMax:
            NewPitch = CurrentPitch + dPitchMax
        elif CurrentPitch > dPitchMax:
            NewPitch = CurrentPitch - dPitchMax
        else:
            NewPitch = 0

    return NewPitch


def NavAngleToMathAngle(NavAngle):
    """
    Pretvara kut iz navigacijskog koordinatnog sustava (heading, track, course) u matematički polarni kut (od X-osi, CCW, radijani).
    Ulaz:
        NavAngle : kut u stupnjevima (može biti bilo koji realan broj)
    Izlaz:
        MathAngle : kut u radijanima (0 = X-osa, CCW)
    """
    # Normalizacija na [0, 360)
    NavAngle = NavAngle - 360 * round(NavAngle / 360) if abs(NavAngle) > 360 else NavAngle
    if NavAngle < 0:
        NavAngle += 360

    # Pretvorba iz navigacijskog u matematički kut
    if NavAngle <= 90:
        NavAngle = 90 - NavAngle
    else:
        NavAngle = 450 - NavAngle

    MathAngle = NavAngle * np.pi / 180  # Pretvori u radijane
    return MathAngle

def azimuth(lat1, lon1, lat2, lon2):
    """
    Vraća azimut (bearing) od točke (lat1, lon1) do (lat2, lon2), u stupnjevima.
    Kompatibilno s MATLAB my_azimuth funkcijom.
    """
    # Pretvori u radijane
    lat1 = lat1 * np.pi / 180
    lon1 = lon1 * np.pi / 180
    lat2 = lat2 * np.pi / 180
    lon2 = lon2 * np.pi / 180
    
    # Izračunaj azimut (ista formula kao MATLAB)
    dlon = lon2 - lon1
    y = np.sin(dlon) * np.cos(lat2)
    x = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)
    az = np.arctan2(y, x)
    
    # Pretvori u stupnjeve i normaliziraj na [0, 360)
    az = (az * 180 / np.pi) % 360
    return az


# Ako želiš azimut u stupnjevima:
# np.degrees(azimuth(lat1, lon1, lat2, lon2))

def crosstrackerror(lat1, lon1, lat2, lon2, latAC, lonAC):
    """
    Izračunava cross-track error (CTE) u nautičkim miljama između točke AC i linije definirane točkama (lat1, lon1) i (lat2, lon2) na sferi.
    Ulazi su u stupnjevima.
    Povrat: cross-track error (najkraća udaljenost od točke do velike kružnice kroz AB), u nautičkim miljama.
    """
    # Pretvori stupnjeve u radijane
    a = np.radians([lat1, lon1])
    b = np.radians([lat2, lon2])
    c = np.radians([latAC, lonAC])

    # Pretvorba iz sfernih u kartezijanske koordinate (MATLAB sph2cart: sph2cart(lon, lat, r))
    def sph2cart(lon, lat, r=1):
        x = r * np.cos(lat) * np.cos(lon)
        y = r * np.cos(lat) * np.sin(lon)
        z = r * np.sin(lat)
        return np.array([x, y, z])

    A = sph2cart(a[1], a[0])
    B = sph2cart(b[1], b[0])
    C = sph2cart(c[1], c[0])

    # Jedinični normalni vektor na ravninu definiranu s A, B i centrom sfere
    n = np.cross(A, B)
    n = n / np.linalg.norm(n)

    # Izračun sinusa kuta theta između vektora OC i ravnine P
    sinTheta = abs(np.dot(C, n))
    theta = np.arcsin(sinTheta)  # u radijanima

    # Pretpostavljeni radijus Zemlje u nautičkim miljama
    radius_nm = 3440.065

    # Povratak površinske udaljenosti (CTE) u nautičkim miljama
    CTE = theta * radius_nm
    return CTE

def matlab_round(n):
    return np.sign(n) * np.floor(np.abs(n) + 0.5)

def BankSet(GP, ACstate, ACcontrol, waypoints, WPTi, ConfMode, Wind, const):
    """
    Bank angle controller. Vraća novi bank angle u radijanima.
    Python verzija usklađena s MATLAB logikom.
    """
    CurrentX = ACstate[0]
    CurrentY = ACstate[1]
    CurrentSpd = ACstate[3]
    CurrentHdg = ACstate[4]
    CurrentBank = ACcontrol[1]

    # Track - azimut iz trenutne pozicije prema idućem waypointu
    # Pretpostavka je da su funkcije my_azimuth i NavAngleToMathAngle prenesene u Python
    DirectTo = azimuth(CurrentY, CurrentX, waypoints[WPTi]['y'], waypoints[WPTi]['x'])
    DirectTo = NavAngleToMathAngle(DirectTo)
    
    Track = DirectTo

    """
    ovo trosi previse vremena, a nije potrebno bar po kodu
    # Cross-track error
    CTE = crosstrackerror(
        waypoints[WPTi-1]['y'], waypoints[WPTi-1]['x'],
        waypoints[WPTi]['y'], waypoints[WPTi]['x'],
        CurrentY, CurrentX
    )
    """

    # Maksimalni dozvoljeni bank angle
    if ConfMode == 'TO':
        BankMax = GP['phi_nom_to'] * np.pi / 180
    else:
        BankMax = GP['phi_nom_oth'] * np.pi / 180

    # Proračuni za presretanje
    WindDir = np.arctan2(Wind[1], Wind[0])
    WindSpd = np.sqrt(Wind[0]**2 + Wind[1]**2)
    """
    WCA = np.arcsin(np.sin(Track - WindDir) * WindSpd / CurrentSpd)
    
    InterceptHdg = DirectTo + WCA

    # Proračun LeadHdg (potpuno odgovara MATLAB logici)
    NSteps = abs(round(CurrentBank / 0.034906585)) """
    asin_arg = np.sin(Track - WindDir) * WindSpd / CurrentSpd
    asin_arg_clipped = np.clip(asin_arg, -1.0, 1.0)
    WCA = np.arcsin(asin_arg_clipped)
    
    InterceptHdg = DirectTo + WCA

    # Korištenje ISPRAVLJENE round funkcije
    NSteps = matlab_round(abs(CurrentBank / 0.034906585))

    sum_turn = 0
    if NSteps > 0:
        if CurrentBank > 0:
            for i in range(int(NSteps) + 1):
                sum_turn += np.tan(CurrentBank - 0.034906585 * i)
        else:
            for i in range(int(NSteps) + 1):
                sum_turn += np.tan(CurrentBank + 0.034906585 * i)
    
    LeadHdg = (const['g0'] * sum_turn / CurrentSpd) if NSteps > 0 else 0

    # Razlika u kursu
    HdgDiff = InterceptHdg - CurrentHdg
    if abs(HdgDiff) > np.pi * 2:
        HdgDiff -= np.pi * 2 * round(HdgDiff / (np.pi * 2))

    DirectionLeft = True if np.sin(HdgDiff) > 0 else False

    # Glavna logika za određivanje novog bank kuta
    NewBank = 0  # Inicijalizacija
    if DirectionLeft and CurrentBank < 0:
        NewBank = CurrentBank + 0.034906585
    elif not DirectionLeft and CurrentBank > 0:
        NewBank = CurrentBank - 0.034906585
    else:
        if DirectionLeft:
            if HdgDiff > LeadHdg:
                NewBank = CurrentBank + 0.034906585
            else:  # HdgDiff <= LeadHdg
                if HdgDiff < 0:
                    NewBank = CurrentBank + 0.034906585
                else:
                    NewBank = CurrentBank - 0.034906585
        else:  # Not DirectionLeft
            if HdgDiff < LeadHdg or HdgDiff > np.pi:
                NewBank = CurrentBank - 0.034906585
            else:  # HdgDiff >= LeadHdg
                NewBank = CurrentBank + 0.034906585

    # ISPRAVAK: Ovaj blok je pomaknut izvan 'else' petlje da odgovara MATLAB logici
    # Ograničavanje na maksimalni bank angle
    if abs(NewBank) > BankMax:
        if NewBank > 0:
            NewBank = NewBank - 0.034906585
        else:
            NewBank = NewBank + 0.034906585
            
    return NewBank


def FuelConsumption(TAS, Thrust, Hp, Cf1, Cf2, Cf3, Cf4, Cfcr, ACtype, PhaseOfFlight):
    """
    Računa potrošnju goriva u kg/s prema fazi leta i tipu motora (BADA 3.15).
    """
    # Pretvorbe jedinica
    TAS = TAS * 1.943844492440605  # m/s u čvorove
    Hp = Hp / 0.3048               # m u ft

    actype = ACtype.lower()
    phase = PhaseOfFlight.lower()

    if actype in ['j', 'jet']:
        eta = Cf1 * (1 + TAS / Cf2)
        Fnom = eta * Thrust / 60000
        Fmin = Cf3 * (1 - Hp / Cf4) / 60
        Fcr = eta * Thrust * Cfcr / 60000
        piston = False
    elif actype in ['t', 'turboprop']:
        eta = Cf1 * (1 - TAS / Cf2) * (TAS / 1000)
        Fnom = eta * Thrust / 60000
        Fmin = Cf3 * (1 - Hp / Cf4) / 60
        Fcr = eta * Thrust * Cfcr / 600000
        piston = False
    elif actype in ['p', 'piston']:
        Fnom = Cf1 / 60
        Fmin = Cf3 / 60
        Fcr = Cf1 * Cfcr / 60
        piston = True
    else:
        Fnom = Fmin = Fcr = -999999999
        piston = False

    if phase in ['climb', 'c']:
        FuelFlow = Fnom
    elif phase in ['cruise', 'l', 'level']:
        FuelFlow = Fcr
    elif phase in ['descent', 'd']:
        if piston:
            FuelFlow = Fmin
        else:
            FuelFlow = max(Fnom, Fmin)
    else:
        FuelFlow = -999999999

    return FuelFlow

def MathtoNavAngle(MathAngle):
    """
    Convert angles such as headings and tracks from math coord system to air nav coord system.
    MathAngle - Angle with x-axis. [radians]
    Air nav angle measured from north clockwise. [degrees]
    """
    NavAngle = MathAngle * 57.29577951308233  # Koristi istu konstantu kao MATLAB
    
    NavAngle = -NavAngle + 90
    
    if abs(NavAngle) > 360:
        NavAngle = NavAngle - 360 * round(NavAngle / 360)  # bounding the value of HdgDiff to +/-360
    
    if NavAngle < 0:
        NavAngle = NavAngle + 360
    
    return NavAngle

#print(MathtoNavAngle(1.569692))
def cart2compass(u, v):
    """
    Pretvara kart. komponente (u,v) u smjer (degN) i brzinu.
    """
    theta, rho = np.arctan2(v, u), np.hypot(u, v)
    theta = np.degrees(theta)
    if theta < 0:
        theta += 360
    # Pretvori u "kompas" kut (0° = sjever, CW)
    if 0 <= theta < 90:
        theta_comp = abs(theta - 90)
    else:
        theta_comp = abs(450 - theta)
    return theta_comp, rho

def ACdrift(course, airspeed, windfrom, windspeed):
    """
    Računa track, groundspeed i wind drift angle.
    course, windfrom: u stupnjevima (navigacijski kutevi, 0° = sjever, CW)
    airspeed, windspeed: u m/s
    """
    windang = np.radians(windfrom - course)
    speedratio = (windspeed / airspeed) * np.sin(windang)
    speedratioC = (windspeed / airspeed) * np.cos(windang)
    if speedratio >= 1 or speedratioC >= 1:
        # Drift correction not possible
        return np.nan, np.nan

    # Track
    track = course - np.degrees(np.arcsin(speedratio))
    track = track % 360  # 0-360
    # Groundspeed
    groundspeed = airspeed * np.sqrt(1 - speedratio ** 2) - windspeed * np.cos(windang)
    # Wind drift angle
    winddriftangle = (track - course + 180) % 360 - 180  # npi2pi
    return track, groundspeed

#HDG=309.27°, ACstate(4)=100.03 m/s, WindDir=90.00°, WindSpd=0.00 m/s
#Track: 309.3°, GS: 100.0 m/s
#print(ACdrift(0.06, 309.39, 90.00, 0.00))
def km2deg(km, radius=6371.0):
    """
    Pretvara udaljenost u kilometrima u stupnjeve luka na sferi (default: Zemlja).
    Parametri:
        km     : udaljenost u kilometrima (može biti float ili array)
        radius : radijus sfere u km (default: 6371.0)
    Povrat:
        deg    : kutna udaljenost u stupnjevima
    """
    return (km / (2 * np.pi * radius)) * 360

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

"""
ACState=2.456195
HDG=309.27°, ACstate(4)=100.03 m/s, WindDir=90.00°, WindSpd=0.00 m/s
Track: 309.3°, GS: 100.0 m/s
"""
"""
print(MathtoNavAngle(1.5707963267948966))
track, gs = ACdrift(309.27, 100.03 , 90.00,0.00)
print(f"track: {track} GS: {gs}") """
def ACStateUpdate(ACstate, ACcontrol, Wind, S, rho, FF, CL, const, method):
    """
    This function updated state of aircraft (ACState) based on aircraft
    control parameters (ACControl), atmosphere data (wind, rho) and aircraft
    parameters.
    """
    ACnewstate = [0] * 6  # Inicijaliziraj kao u MATLAB-u

    # converting wind components to Wind direction (WindDir) and Wind speed (WindSpd)
    WindDir, WindSpd = cart2compass(-Wind[0], -Wind[1])
    # windDir must be reversed because this functions just sum v and u vector
    # that why components are multiplied by -1

    # calculate track and groundspeed (GS) using wind and AC state data
    #print(f"\nACState={ACstate[4]:.6f}")
    HDG = MathtoNavAngle(ACstate[4])  # ACstate(5) u MATLAB-u je ACstate[4] u Python-u
    #print(f'HDG={HDG:.2f}°, ACstate(3)={ACstate[3]:.2f} m/s, WindDir={WindDir:.2f}°, WindSpd={WindSpd:.2f} m/s')

    track, GS = ACdrift(HDG, ACstate[3], WindDir, WindSpd)  # ACstate(4) u MATLAB-u je ACstate[3] u Python-u
    #print(f'Track: {track:.1f}°, GS: {GS:.1f} m/s')
    
    # GS is in (m/s) and step of simulation is 1 second so it is required to
    # change it to km/s for km2deg to work.
    GS = GS / 1000  # m/s -> km/s
    dist = km2deg(GS)  # converting distance covered in 1 second to degrees

    # reckon is function to calculate new position using current position,
    # distance and direction of flight. Result is [lat,lon] of new point
    newpoint = reckon(ACstate[1], ACstate[0], dist, track)  # ACstate(2), ACstate(1) u MATLAB-u

    # update AC x position
    ACnewstate[0] = newpoint[1]  # ACnewstate(1) = newpoint(2) u MATLAB-u
    # update AC y position  
    ACnewstate[1] = newpoint[0]  # ACnewstate(2) = newpoint(1) u MATLAB-u
    # update AC h position
    ACnewstate[2] = ACstate[2] + ACstate[3] * np.sin(ACcontrol[2]) + Wind[2]  # ACstate(3)+ACstate(4)*sin(ACcontrol(3))+Wind(3)
    # update AC TAS
    
    """
    if method == "tg":
        print(f'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA    ACstate[3]={ACstate[3]:}, ACstate[5]={ACstate[5]:}, ACcontrol[0]={ACcontrol[0]:}, ACcontrol[2]={ACcontrol[2]:}, ACcontrol[3]={ACcontrol[3]:}, S={S:}, rho={rho:}, g0={const["g0"]:}')
        print(f'((ACcontrol[3] * S * rho * ACstate[3]**2) / (2 * ACstate[5]))={((ACcontrol[3] * S * rho * (ACstate[3]**2)) / (2 * ACstate[5]))}')"""

    ACnewstate[3] = (ACstate[3] - 
                     ((ACcontrol[3] * S * rho * ACstate[3]**2) / (2 * ACstate[5])) - 
                     const['g0'] * np.sin(ACcontrol[2]) + 
                     ACcontrol[0] / ACstate[5])  # ACstate(4), ACcontrol(4), ACcontrol(3), ACcontrol(1), ACstate(6) u MATLAB-u
    # update AC hdg
    ACnewstate[4] = ACstate[4] + (CL * S * rho * ACstate[3] * np.sin(ACcontrol[1])) / (2 * ACstate[5])  # ACstate(5), ACstate(4), ACcontrol(2), ACstate(6) u MATLAB-u
    #print(f"Heading update: ACstate[4]={ACstate[4]:.6f}, CL={CL:.6f}, S={S:.2f}, rho={rho:.6f}, TAS={ACstate[3]:.2f}, bank={ACcontrol[1]:.6f}, mass={ACstate[5]:.2f}, delta={(CL * S * rho * ACstate[3] * np.sin(ACcontrol[1])) / (2 * ACstate[5]):.6f}")

    # update AC mass
    ACnewstate[5] = ACstate[5] - FF  # ACstate(6) u MATLAB-u

    #print(f'newpoint: [{newpoint[0]:.6f}, {newpoint[1]:.6f}]')
    #print(f'ACnewstate: [{ACnewstate[0]:.4f} {ACnewstate[1]:.4f} {ACnewstate[2]:.4f} {ACnewstate[3]:.4f} {ACnewstate[4]:.4f} {ACnewstate[5]:.4f}]')
    #print()
    return ACnewstate


#print(1 + 70 * np.sin(1) + 0)  # ACstate(3)+ACstate(4)*sin(ACcontrol(3))+Wind(3))

def CDA(waypoints, ACstate, WPTi, dT, const, ACmode, GP, apfdata, opsdata, ACcontrol, Wind):
    """
    this function will check if AC is landing at next WPT and calculate TOD
    for descent
    """
    
    t = 1  # MATLAB počinje s t=1, ne t=0
    D = []
    
    # Postavi početno stanje na maksimalnu visinu rute
    ACstate[2] = max(wp['z'] for wp in waypoints)  # ACstate(3) u MATLAB-u
    ACmode['CL'] = 'C'
    ACmode['ACC'] = 'C'
    ACmode['ConfigMode'] = 'CL'
    ACstate[3] = opsdata['Vdesref']  # ACstate(4) u MATLAB-u
    #print(f"UNTUAR CDA 4: [{ACstate[3]}]")
    ACstate[5] = ACstate[5] * 0.9  # ACstate(6) u MATLAB-u
    #print(f"UNTUAR CDA 6: [{ACstate[5]}]")

    # Prvi silazak
    while ACstate[2] > waypoints[-1]['z'] + 150:  # waypoints(end).z u MATLAB-u
        
        # 1. Determine conditions at current alt for current meteo
        AtmHp_T, AtmHp_p, AtmHp_rho, AtmHp_a = AtmosphereAtHp(ACstate[2], dT, const)
        AtmHp = {'T': AtmHp_T, 'p': AtmHp_p, 'rho': AtmHp_rho, 'a': AtmHp_a}

        # 2. Determine altitude
        if ACstate[2] > 11000:
            ACmode['tropo'] = True
        else:
            ACmode['tropo'] = False

        # 3. Determine all speeds (TAS to CAS to Mach)
        CAScurrent = TAStoCAS(AtmHp['p'], AtmHp['rho'], ACstate[3])  # ACstate(4) - TAS u MATLAB-u
        MACHcurrent = TAStoMach(AtmHp['T'], ACstate[3])
        TransAlt = CrossAlt(CAScurrent, MACHcurrent, dT, const)
        
        if ACstate[2] > TransAlt:  # switch for constant CAS or Mach in respect to Transition (crossover) altitude
            ACmode['SpeedMode'] = 'M'
        else:
            ACmode['SpeedMode'] = 'C'

        # 4. Determine phase of flight: (T)ake-(O)ff = TO, (I)nitial (C)limb, (CL)ean, (APP)roach, (L)an(D)in(G)
        #print(f"CDAAAAAAAAAAAAAAAAAAAAAA: ACstate[2]={ACstate[2]}, waypoints[WPTi]['z']={waypoints[WPTi]['z']}, ACmode['CL']={ACmode['CL']}, GP['Htol']={GP['Htol']}",flush=True)
        ACmode['CL'] = CLModeSet(ACstate[2], waypoints[-1]['z'], ACmode['CL'], GP['Htol'])  # setting AC regime
        #print(f"CDAAAAAAAAAAAAAAAAAAAAAA: ACstate[2]={ACstate[2]}, CAScurrent={CAScurrent}, ACmode['CL']={ACmode['CL']}",flush=True)
        ACmode['ConfigMode'] = ConfigModeSet(GP, opsdata, ACstate[2], CAScurrent, ACmode['CL'])  # Setting config mode
        #print(f"CDAAAAAAAAAAAAAAAAAAAAAA: ACmode['ConfigMode'] = {ACmode['ConfigMode']}",flush=True)
        #print()
        # 5. Determine acceleration mode: (A)ccelerating, (C)onstant, (D)ecelerating, and desired TAS
        desiredTAS = DesiredTASf(GP, apfdata, opsdata, ACmode['CL'], ACstate[5], AtmHp, ACmode['ConfigMode'], ACstate[2], TransAlt, const)
        ACmode['ACC'] = ACCModeSet(ACstate[3], desiredTAS, ACmode['ACC'], GP['Vtol'])

        # 6. Determine energy share factor
        esf = energysf(MACHcurrent, AtmHp['T'], dT, ACmode['ACC'], ACmode['CL'], ACmode['SpeedMode'], ACmode['tropo'], const)

        # 7. Determine ACcontrol inputs
        #   7.a. Lift and Drag
        CL = LiftCoefficient(ACstate[5], AtmHp['rho'], opsdata['wingsurf'], ACstate[3], ACcontrol[1], const)  # ACcontrol(2) u MATLAB-u
        ACcontrol[3] = ACDragCoef(opsdata, ACmode['ConfigMode'], CL)  # ACcontrol(4) u MATLAB-u
        #   7.b. Thrust
        ACcontrol[0] = ThrustSet(GP, opsdata, ACstate, ACcontrol, AtmHp, ACmode, desiredTAS, const)  # ACcontrol(1) u MATLAB-u
        
        #   7.c. Pitch
        #print(f"CDA::::::::::::::::::::::::ACcontrol prije poziva pitchset = {ACcontrol[2]}")
        ACcontrol[2] = PitchSet(opsdata, GP, ACstate, ACcontrol, waypoints[-1]['z'], esf, AtmHp, const)  # ACcontrol(3) u MATLAB-u
        
        # print(f"CDA::::::::::::::::::::::::ACcontrol polsije poziva pitchset = {ACcontrol[2]}")
        #   7.d. Bank
        ACcontrol[1] = BankSet(GP, ACstate, ACcontrol, waypoints, WPTi, ACmode['ConfigMode'], Wind, const)  # ACcontrol(2) u MATLAB-u

        # 8. Fuel consumption
        FF = FuelConsumption(ACstate[3], ACcontrol[0], ACstate[2],  # ACstate(4), ACcontrol(1), ACstate(3) u MATLAB-u
            opsdata['Cf1'], opsdata['Cf2'], opsdata['Cf3'], opsdata['Cf4'],
            opsdata['Cfcr'], opsdata['engtype'], ACmode['CL'])

        # calculate next ACstate
        ACSm1 = ACstate.copy()
        ACstate = ACStateUpdate(ACstate, ACcontrol, Wind, opsdata['wingsurf'], AtmHp['rho'], FF, CL, const,"cda")
        
        D.append(distance(ACstate[1], ACstate[0], ACSm1[1], ACSm1[0]))  # ZADRŽANO distance
        t = t + 1

    Ds = sum(D)
    p1 = len(waypoints)  # size(waypoints,2) u MATLAB-u
    leg = distance(waypoints[p1-1]['y'], waypoints[p1-1]['x'], waypoints[p1-2]['y'], waypoints[p1-2]['x'])  # ZADRŽANO distance
    maxleg = 0
    for m in range(1, p1):  # for m=2:p1 u MATLAB-u
        maxleg = maxleg + distance(waypoints[m]['y'], waypoints[m]['x'], waypoints[m-1]['y'], waypoints[m-1]['x'])  # ZADRŽANO distance

    if maxleg > Ds:  # if ac has passed its CDA mark before start of simulation
        while Ds > leg:
            Ds = Ds - leg
            p1 = p1 - 1
            leg = distance(waypoints[p1-1]['y'], waypoints[p1-1]['x'], waypoints[p1-2]['y'], waypoints[p1-2]['x'])  # ZADRŽANO distance

        az = azimuth(waypoints[p1-2]['y'], waypoints[p1-2]['x'], waypoints[p1-1]['y'], waypoints[p1-1]['x'])  # ZADRŽANO azimuth
        az = az - 180
        if az < 0:
            az = az + 360
    else:
        waypoints[0]['name'] = 'CDA'

    # if waypoints(p1-1).z < max([waypoints(:).z]) u MATLAB-u
    if waypoints[p1-2]['z'] < max(wp['z'] for wp in waypoints):
        #print("TU SAM")
        
        ACstate[2] = waypoints[p1-2]['z']  # ACstate(3) u MATLAB-u
        
        while ACstate[2] > waypoints[-1]['z'] + 150:

            # 1. Determine conditions at current alt for current meteo
            AtmHp_T, AtmHp_p, AtmHp_rho, AtmHp_a = AtmosphereAtHp(ACstate[2], dT, const)
            AtmHp = {'T': AtmHp_T, 'p': AtmHp_p, 'rho': AtmHp_rho, 'a': AtmHp_a}

            # 2. Determine altitude
            if ACstate[2] > 11000:
                ACmode['tropo'] = True
            else:
                ACmode['tropo'] = False
            
            # 3. Determine all speeds (TAS to CAS to Mach)
            CAScurrent = TAStoCAS(AtmHp['p'], AtmHp['rho'], ACstate[3])  # ACstate(4) - TAS
            MACHcurrent = TAStoMach(AtmHp['T'], ACstate[3])
            TransAlt = CrossAlt(CAScurrent, MACHcurrent, dT, const)

            if ACstate[2] > TransAlt:  # switch for constant CAS or Mach in respect to Transition (crossover) altitude
                ACmode['SpeedMode'] = 'M'
            else:
                ACmode['SpeedMode'] = 'C'

            # 4. Determine phase of flight: (T)ake-(O)ff = TO, (I)nitial (C)limb, (CL)ean, (APP)roach, (L)an(D)in(G)
            #print(f"CDAAAAAAAAAAAAAAAAAAAAAA: ACstate[2]={ACstate[2]}, waypoints[WPTi]['z']={waypoints[WPTi]['z']}, ACmode['CL']={ACmode['CL']}, GP['Htol']={GP['Htol']}",flush=True)
            ACmode['CL'] = CLModeSet(ACstate[2], waypoints[-1]['z'], ACmode['CL'], GP['Htol'])  # setting AC regime
            #print(f"CDAAAAAAAAAAAAAAAAAAAAAA: ACstate[2]={ACstate[2]}, CAScurrent={CAScurrent}, ACmode['CL']={ACmode['CL']}",flush=True)
            ACmode['ConfigMode'] = ConfigModeSet(GP, opsdata, ACstate[2], CAScurrent, ACmode['CL'])  # Setting config mode
            #print(f"CDAAAAAAAAAAAAAAAAAAAAAA: ACmode['ConfigMode'] = {ACmode['ConfigMode']}",flush=True)
            #print()

            # 5. Determine acceleration mode: (A)ccelerating, (C)onstant, (D)ecelerating, and desired TAS
            desiredTAS = DesiredTASf(GP, apfdata, opsdata, ACmode['CL'], ACstate[5], AtmHp, ACmode['ConfigMode'], ACstate[2], TransAlt, const)
            ACmode['ACC'] = ACCModeSet(ACstate[3], desiredTAS, ACmode['ACC'], GP['Vtol'])

            # 6. Determine energy share factor
            esf = energysf(MACHcurrent, AtmHp['T'], dT, ACmode['ACC'], ACmode['CL'], ACmode['SpeedMode'], ACmode['tropo'], const)

            # 7. Determine ACcontrol inputs
            #   7.a. Lift and Drag
            CL = LiftCoefficient(ACstate[5], AtmHp['rho'], opsdata['wingsurf'], ACstate[3], ACcontrol[1], const)
            ACcontrol[3] = ACDragCoef(opsdata, ACmode['ConfigMode'], CL)
            #   7.b. Thrust
            
            ACcontrol[0] = ThrustSet(GP, opsdata, ACstate, ACcontrol, AtmHp, ACmode, desiredTAS, const)
            
            #   7.c. Pitch
            #print(f"CDA drugi pozivvvvvvvv::::::::::::::::::::::::ACcontrol prije poziva pitchset = {ACcontrol[2]}")
            ACcontrol[2] = PitchSet(opsdata, GP, ACstate, ACcontrol, waypoints[-1]['z'], esf, AtmHp, const)
            #print(f"CDA drugi pozivvvvvvvv::::::::::::::::::::::::ACcontrol polsije poziva pitchset = {ACcontrol[2]}")
            #   7.d. Bank
            ACcontrol[1] = BankSet(GP, ACstate, ACcontrol, waypoints, WPTi, ACmode['ConfigMode'], Wind, const)

            # 8. Fuel consumption
            FF = FuelConsumption(ACstate[3], ACcontrol[0], ACstate[2],
                opsdata['Cf1'], opsdata['Cf2'], opsdata['Cf3'], opsdata['Cf4'],
                opsdata['Cfcr'], opsdata['engtype'], ACmode['CL'])

            # calculate next ACstate
            ACSm1 = ACstate.copy()
            ACstate = ACStateUpdate(ACstate, ACcontrol, Wind, opsdata['wingsurf'], AtmHp['rho'], FF, CL, const,"cda")

            D.append(distance(ACstate[1], ACstate[0], ACSm1[1], ACSm1[0]))  # ZADRŽANO distance
            t = t + 1
        
        Ds = sum(D)
        p1 = len(waypoints)
        leg = distance(waypoints[p1-1]['y'], waypoints[p1-1]['x'], waypoints[p1-2]['y'], waypoints[p1-2]['x'])  # ZADRŽANO distance
        
        for m in range(1, p1):  # for m=2:p1 u MATLAB-u
            maxleg = maxleg + distance(waypoints[m]['y'], waypoints[m]['x'], waypoints[m-1]['y'], waypoints[m-1]['x'])  # ZADRŽANO distance
    
        if maxleg > Ds:

            while Ds > leg:
                Ds = Ds - leg
                p1 = p1 - 1
                if p1 == 1:  # if p1==1, p1=2; u MATLAB-u
                    p1 = 2
                leg = distance(waypoints[p1-1]['y'], waypoints[p1-1]['x'], waypoints[p1-2]['y'], waypoints[p1-2]['x'])  # ZADRŽANO distance

            az = azimuth(waypoints[p1-2]['y'], waypoints[p1-2]['x'], waypoints[p1-1]['y'], waypoints[p1-1]['x'])  # ZADRŽANO azimuth
            az = az - 180
            if az < 0:
                az = az + 360
        else:
            waypoints[0]['name'] = 'CDA'

    if waypoints[0]['name'] != 'CDA':  # if ~strcmp(waypoints(1).name,'CDA') u MATLAB-u
        latout, lonout = reckon(waypoints[p1-1]['y'], waypoints[p1-1]['x'], Ds, az)  # ZADRŽANO reckon

        A = waypoints[:p1-1]  # waypoints(1:p1-1) u MATLAB-u
        B = {
            'y': latout,
            'x': lonout,
            'z': waypoints[p1-2]['z'],  # waypoints(p1-1).z u MATLAB-u
            'flyover': 0,
            'hist': 1,
            'name': 'CDA'
        }
        C = waypoints[p1-1:]  # waypoints(p1:end) u MATLAB-u

        wptlist = A + [B] + C  # [A,B,C] u MATLAB-u
        # [wptlist(p1+1:end).z]=deal(0); set all alts to 0 after CDA mark
        for i in range(p1, len(wptlist)):
            wptlist[i]['z'] = 0
    else:
        wptlist = waypoints
    import json
    #print("KAD ZAVRSI U CDA")
    #print(json.dumps(ACstate, indent=4, default=str))
    return wptlist

