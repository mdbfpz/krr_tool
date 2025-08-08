import numpy as np
const = {
    'earthradius': 6325766,
    'Hp_trop': 11000,
    'T0': 288.15,
    'P0': 101325,
    'g0': 9.8067,
    'a0': 340.294,
    'Tisa_trop': 216.65,
    'BetaT': -0.0065,
    'kappa': 1.4,
    'R': 287.0529,
    'ft': 0.3048
}

def AtmosphereAtHp(Hp, dT, const):
    """
    Izračun atmosferskih uvjeta na zadanoj geopotencijalnoj visini.
    
    Ulazi:
        Hp    - geopotencijalna visina (u metrima)
        dT    - temperaturna razlika od ISA standarda na MSL (u Kelvinima)
        const - rječnik ili objekt s konstantama:
                T0, P0, g0, BetaT, Hp_trop, R, Tisa_trop, kappa
    
    Izlazi:
        T   - temperatura na visini (K)
        p   - tlak (Pa)
        rho - gustoća zraka (kg/m^3)
        a   - brzina zvuka (m/s)
    """
    # Troposfera
    if Hp < const['Hp_trop']:
        T = const['T0'] + dT + const['BetaT'] * Hp
        p = const['P0'] * ((T - dT) / const['T0']) ** (-const['g0'] / (const['BetaT'] * const['R']))
    # Iznad tropopauze
    else:
        T = const['T0'] + dT + const['BetaT'] * const['Hp_trop']
        PatHptrop = const['P0'] * (1 + (const['BetaT'] * const['Hp_trop']) / const['T0']) ** (-const['g0'] / (const['BetaT'] * const['R']))
        p = PatHptrop * np.exp((-const['g0'] / (const['R'] * const['Tisa_trop'])) * (Hp - const['Hp_trop']))
    
    rho = AirDensity(const, p, T)
    a = np.sqrt(const['kappa'] * const['R'] * T)
    return T, p, rho, a

def AirDensity(const, p, T):
    """
    Gustoća zraka prema idealnom plinskom zakonu.
    """
    return p / (const['R'] * T)
