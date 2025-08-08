def read_APF_file(apffilename):
    """
    Parsira BADA APF datoteku i vraća strukturu s podacima.
    
    Ulaz:
        apffilename - naziv APF datoteke koju treba parsirati
    
    Izlaz:
        apfdata - rječnik sa svim podacima pročitanim iz APF datoteke
    """
    
    # Inicijalizacija rječnika za pohranu podataka
    apfdata = {}
    
    # Otvaranje datoteke za čitanje
    with open(apffilename, 'r') as fid:
        # Čitanje svih linija datoteke
        data = fid.readlines()
        
        # Parsiranje LO i HI vrijednosti iz 13. linije
        a = data[12]  # Indeksiranje počinje od 0 u Pythonu
        apfdata['LO'] = a[9:15].strip()
        apfdata['HI'] = float(a[68:74].strip())
        
        # Inicijalizacija ugniježđenih rječnika
        apfdata['V_cl1'] = {}
        apfdata['V_cl2'] = {}
        apfdata['M_cl'] = {}
        apfdata['V_cr1'] = {}
        apfdata['V_cr2'] = {}
        apfdata['M_cr'] = {}
        apfdata['M_des'] = {}
        apfdata['V_des2'] = {}
        apfdata['V_des1'] = {}
        
        # Parsiranje LO vrijednosti iz 21. linije
        b_line = data[20]
        b_values = b_line[27:].split()
        apfdata['V_cl1']['LO'] = float(b_values[0])
        apfdata['V_cl2']['LO'] = float(b_values[1])
        apfdata['M_cl']['LO'] = float(b_values[2])
        apfdata['V_cr1']['LO'] = float(b_values[3])
        apfdata['V_cr2']['LO'] = float(b_values[4])
        apfdata['M_cr']['LO'] = float(b_values[5])
        apfdata['M_des']['LO'] = float(b_values[6])
        apfdata['V_des2']['LO'] = float(b_values[7])
        apfdata['V_des1']['LO'] = float(b_values[8])
        
        # Parsiranje AV vrijednosti iz 22. linije
        c_line = data[21]
        c_values = c_line[27:].split()
        apfdata['V_cl1']['AV'] = float(c_values[0])
        apfdata['V_cl2']['AV'] = float(c_values[1])
        apfdata['M_cl']['AV'] = float(c_values[2])
        apfdata['V_cr1']['AV'] = float(c_values[3])
        apfdata['V_cr2']['AV'] = float(c_values[4])
        apfdata['M_cr']['AV'] = float(c_values[5])
        apfdata['M_des']['AV'] = float(c_values[6])
        apfdata['V_des2']['AV'] = float(c_values[7])
        apfdata['V_des1']['AV'] = float(c_values[8])
        
        # Parsiranje HI vrijednosti iz 23. linije
        d_line = data[22]
        d_values = d_line[27:].split()
        apfdata['V_cl1']['HI'] = float(d_values[0])
        apfdata['V_cl2']['HI'] = float(d_values[1])
        apfdata['M_cl']['HI'] = float(d_values[2])
        apfdata['V_cr1']['HI'] = float(d_values[3])
        apfdata['V_cr2']['HI'] = float(d_values[4])
        apfdata['M_cr']['HI'] = float(d_values[5])
        apfdata['M_des']['HI'] = float(d_values[6])
        apfdata['V_des2']['HI'] = float(d_values[7])
        apfdata['V_des1']['HI'] = float(d_values[8])
    
    return apfdata

apfdata = read_APF_file("test.APF")
#print(apfdata)

def read_OPF_file(opsfilename):
    """
    Parsira BADA OPF datoteku i vraća rječnik s podacima.
    
    Ulaz:
        opsfilename - naziv OPF datoteke koju treba parsirati
    
    Izlaz:
        opsdata - rječnik sa svim podacima pročitanim iz OPF datoteke
    """
    
    # Inicijalizacija rječnika za pohranu podataka
    opsdata = {}
    
    # Otvaranje datoteke za čitanje
    with open(opsfilename, 'r') as fid:
        # Čitanje svih linija datoteke
        data = fid.readlines()
        
        # Parsiranje tipa zrakoplova, broja motora, tipa motora i kategorije turbulencije
        q_line = data[13]  # Indeks 13 odgovara liniji 14 u MATLAB-u
        q_parts = q_line[5:].split()
        opsdata['actype'] = q_parts[0]
        opsdata['numengines'] = int(q_parts[1])
        opsdata['engtype'] = q_parts[3]
        opsdata['wakecat'] = q_parts[4]
        
        # Parsiranje referentne, minimalne, maksimalne mase, korisnog tereta i težine
        a_line = data[18]  # Indeks 18 odgovara liniji 19 u MATLAB-u
        a_values = a_line[6:].split()
        opsdata['mref'] = float(a_values[0])
        opsdata['mmin'] = float(a_values[1])
        opsdata['mmax'] = float(a_values[2])
        opsdata['mpyld'] = float(a_values[3])
        opsdata['gw'] = float(a_values[4])
        
        # Parsiranje operativnih ograničenja
        b_line = data[21]  # Indeks 21 odgovara liniji 22 u MATLAB-u
        b_values = b_line[6:].split()
        opsdata['VMO'] = float(b_values[0])
        opsdata['MMO'] = float(b_values[1])
        opsdata['maxalt'] = float(b_values[2])
        opsdata['Hmax'] = float(b_values[3])
        opsdata['tempgrad'] = float(b_values[4])
        
        # Parsiranje aerodinamičkih parametara
        c_line = data[25]  # Indeks 25 odgovara liniji 26 u MATLAB-u
        c_values = c_line[6:].split()
        opsdata['wingsurf'] = float(c_values[0])
        opsdata['Clbo'] = float(c_values[1])
        opsdata['k'] = float(c_values[2])
        #print(opsdata['k'])
        opsdata['CM16'] = float(c_values[3])
        
        # Inicijalizacija ugniježđenih rječnika
        opsdata['Vstall'] = {}
        opsdata['Cd0'] = {}
        opsdata['Cd2'] = {}
        
        # Parsiranje aerodinamičkih koeficijenata za različite konfiguracije
        d_line = data[28]  # Indeks 28 odgovara liniji 29 u MATLAB-u
        d_values = d_line[17:].split()
        opsdata['Vstall']['CR'] = float(d_values[0])
        opsdata['Cd0']['CR'] = float(d_values[1])
        opsdata['Cd2']['CR'] = float(d_values[2])
        
        e_line = data[29]  # Indeks 29 odgovara liniji 30 u MATLAB-u
        e_values = e_line[17:].split()
        opsdata['Vstall']['IC'] = float(e_values[0])
        opsdata['Cd0']['IC'] = float(e_values[1])
        opsdata['Cd2']['IC'] = float(e_values[2])
        
        f_line = data[30]  # Indeks 30 odgovara liniji 31 u MATLAB-u
        f_values = f_line[17:].split()
        opsdata['Vstall']['TO'] = float(f_values[0])
        opsdata['Cd0']['TO'] = float(f_values[1])
        opsdata['Cd2']['TO'] = float(f_values[2])
        
        g_line = data[31]  # Indeks 31 odgovara liniji 32 u MATLAB-u
        g_values = g_line[17:].split()
        opsdata['Vstall']['AP'] = float(g_values[0])
        opsdata['Cd0']['AP'] = float(g_values[1])
        opsdata['Cd2']['AP'] = float(g_values[2])
        
        h_line = data[32]  # Indeks 32 odgovara liniji 33 u MATLAB-u
        h_values = h_line[17:].split()
        opsdata['Vstall']['LD'] = float(h_values[0])
        opsdata['Cd0']['LD'] = float(h_values[1])
        opsdata['Cd2']['LD'] = float(h_values[2])
        
        # Parsiranje koeficijenta otpora za spušteni stajni trap
        i_line = data[38]  # Indeks 38 odgovara liniji 39 u MATLAB-u
        i_values = i_line[28:].split()
        opsdata['Cd0']['geardown'] = float(i_values[0])
        
        # Parsiranje koeficijenata potiska
        j_line = data[44]  # Indeks 44 odgovara liniji 45 u MATLAB-u
        j_values = j_line[6:].split()
        opsdata['Ctc1'] = float(j_values[0])
        opsdata['Ctc2'] = float(j_values[1])
        opsdata['Ctc3'] = float(j_values[2])
        opsdata['Ctc4'] = float(j_values[3])
        opsdata['Ctc5'] = float(j_values[4])
        
        # Parsiranje koeficijenata potiska za spuštanje
        k_line = data[46]  # Indeks 46 odgovara liniji 47 u MATLAB-u
        k_values = k_line[6:].split()
        opsdata['Ctdeslow'] = float(k_values[0])
        opsdata['Ctdeshi'] = float(k_values[1])
        opsdata['Hpdes'] = float(k_values[2])
        opsdata['Ctdesapp'] = float(k_values[3])
        opsdata['Ctdesld'] = float(k_values[4])
        
        # Parsiranje referentnih brzina za spuštanje
        l_line = data[48]  # Indeks 48 odgovara liniji 49 u MATLAB-u
        l_values = l_line[6:].split()
        opsdata['Vdesref'] = float(l_values[0])
        opsdata['Mdesref'] = float(l_values[1])
        
        # Parsiranje koeficijenata potrošnje goriva
        m_line = data[51]  # Indeks 51 odgovara liniji 52 u MATLAB-u
        m_values = m_line[6:].split()
        opsdata['Cf1'] = float(m_values[0])
        opsdata['Cf2'] = float(m_values[1])
        
        n_line = data[53]  # Indeks 53 odgovara liniji 54 u MATLAB-u
        n_values = n_line[6:].split()
        opsdata['Cf3'] = float(n_values[0])
        opsdata['Cf4'] = float(n_values[1])
        
        o_line = data[55]  # Indeks 55 odgovara liniji 56 u MATLAB-u
        o_values = o_line[6:].split()
        opsdata['Cfcr'] = float(o_values[0])
        
        # Parsiranje fizičkih dimenzija zrakoplova
        p_line = data[58]  # Indeks 58 odgovara liniji 59 u MATLAB-u
        p_values = p_line[6:].split()
        opsdata['TOL'] = float(p_values[0])
        opsdata['LDL'] = float(p_values[1])
        opsdata['span'] = float(p_values[2])
        opsdata['length'] = float(p_values[3])
    
    return opsdata

opsdata = read_OPF_file("testOPF.OPF")
#print(opsdata)