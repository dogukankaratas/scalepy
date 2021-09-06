import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from scipy import stats
from scipy import integrate
import scipy.signal as signal
from scipy.signal import butter, lfilter
from scipy.interpolate import interp1d


def select_records(magnitude_range, vs_range, rjb, fault_mechanism):


    ### Kayıtların seçilmesi ###

    ## Dosyaların okunması
    FlatFile_path_name = "C:\\Users\\KaratasD\\Desktop\\PRJ_2729\\EQE_Selection_Scaling\\0-Data\\flat_final.xlsx"
    

    eqe_df = pd.read_excel(FlatFile_path_name)

    ## Inputların alınması
    min_m, max_m = [float(x) for x in magnitude_range.split()]
    min_vs, max_vs = [float(x) for x in vs_range.split()]
    min_d, max_d = [float(x) for x in rjb.split()]



    eqe_s = eqe_df[(eqe_df["Earthquake Magnitude"] > min_m) & (eqe_df["Earthquake Magnitude"] < max_m) 
                      & (eqe_df["Vs30 (m/s) selected for analysis"] > min_vs) & (eqe_df["Vs30 (m/s) selected for analysis"] < max_vs)
                      &  (eqe_df["EpiD (km)"] > min_d) & (eqe_df["EpiD (km)"] < max_d)]


    fay = fault_mechanism

    if fay == 'Strike - Slip' :
        eqe_selected = eqe_s[(eqe_s["Rake Angle (deg)"] > -180) & (eqe_s["Rake Angle (deg)"] < -150) | (eqe_s["Rake Angle (deg)"] > -30) & (eqe_s["Rake Angle (deg)"] < 30) | 
                         (eqe_s["Rake Angle (deg)"] > 150) & (eqe_s["Rake Angle (deg)"] < 180)]
    
    elif fay ==  'Normal' :
        eqe_selected = eqe_s[(eqe_s["Rake Angle (deg)"] > -120) & (eqe_s["Rake Angle (deg)"] < -60)]
                         
    elif fay ==  'Reverse' :
        eqe_selected = eqe_s[(eqe_s["Rake Angle (deg)"] > 60) & (eqe_s["Rake Angle (deg)"] < 120)]

    elif fay ==  'Reverse - Oblique' :
        eqe_selected = eqe_s[(eqe_s["Rake Angle (deg)"] > 30) & (eqe_s["Rake Angle (deg)"] < 60) | (eqe_s["Rake Angle (deg)"] > 120) & (eqe_s["Rake Angle (deg)"] < 150)]

    elif fay ==  'Normal - Oblique' :
        eqe_selected = eqe_s[(eqe_s["Rake Angle (deg)"] > -150) & (eqe_s["Rake Angle (deg)"] < -120) | (eqe_s["Rake Angle (deg)"] > -60) & (eqe_s["Rake Angle (deg)"] < -30)]

    else : 
        print("Invalid Mechanism!")
    
    
    ##Seçilen Kayıtların 11'e indirgenmesi

    eqe_records = eqe_selected.sample(n=11)
    eqe_records.drop(eqe_records.columns.difference(['Record Sequence Number', 'Earthquake Name', 
                                                 'File Name (Horizontal 1)', 'File Name (Horizontal 2)']), 1, inplace = True)
    eqe_records = eqe_records.rename(columns={'Record Sequence Number': 'ID'})

    del eqe_selected
    del eqe_s
    del eqe_df

    return eqe_records


def ResponseSpectra(acceleration, sampling_interval, damping_ratio = 0.05, Tmax = 4, T_int = 0.05):

    '''      
    Response spectra using piecewise
    
    Input:
        T: vector with periods (s)
        s: acceleration time series
        zi: damping ratio
        dt: time steps for s
    
    Returns:
        PSA, PSV, SA, SV, SD
    
    '''

    import numpy as np
    T = np.arange(  T_int , Tmax + T_int , T_int)
    s = acceleration
    zi = damping_ratio
    dt = sampling_interval

    pi = np.pi
    
    nper = np.size(T)						      # number of natural periods
    n    = np.size(s)                             # length of record
    
    SD   = np.zeros(nper)				          # rel. displac. spectrum
    SV   = np.zeros(nper)				          # rel. vel. spectrum
    SA   = np.zeros(nper)				          # total acc. spectrum	
     
    
    for k in range(nper):
       wn = 2*pi/T[k]
       wd = wn*(1-zi**2)**(1/2)
       
       u = np.zeros((2,n))          # matrix with velocities and displacements
       
       ex = np.exp(-zi*wn*dt)
       cwd = np.cos(wd*dt)
       swd = np.sin(wd*dt)
       zisq = 1/(np.sqrt(1-(zi**2)))
    
       a11 = ex*(cwd+zi*zisq*swd)
       a12 = (ex/wd)*swd
       a21 = -wn*zisq*ex*swd
       a22 = ex*(cwd-zi*zisq*swd)
    
       b11 = ex*(((2*zi**2-1)/((wn**2)*dt)+zi/wn)*(1/wd)*np.sin(wd*dt)+
           (2*zi/((wn**3)*dt)+1/(wn**2))*np.cos(wd*dt))-2*zi/((wn**3)*dt)
       b12 = -ex*(((2*zi**2-1)/((wn**2)*dt))*(1/wd)*np.sin(wd*dt)+
           (2*zi/((wn**3)*dt))*np.cos(wd*dt))-(1/(wn**2))+2*zi/((wn**3)*dt)
       b21 = -((a11-1)/((wn**2)*dt))-a12
       b22 = -b21-a12
       
       A = np.array([[a11,a12],[a21,a22]])
       B = np.array([[b11,b12],[b21,b22]])
    
       for q in range(n-1):
          u[:,q+1] = np.dot(A,u[:,q]) + np.dot(B,np.array([s[q],s[q+1]]))
       
       at = -2*wn*zi*u[1,:]-(wn**2)*u[0,:]
       
       SD[k]   = np.max( np.abs(u[0,:]) )
       SV[k]   = np.max( np.abs(u[1,:]) )
       SA[k]   = np.max( np.abs(at) )
    
    PSV = (2*pi/T)*SD                    # pseudo-vel. spectrum
    PSA = (2*pi/T)**2 *SD  	             # pseudo-accel. spectrum
    
    return T , SD, SV, SA   


def ReadPeer(filepath):

    fileName = filepath.split(os.path.sep)[-1].replace(".at2", "")
    with open(filepath) as f:
        acc = []
        for no, line in enumerate(f.readlines()):
            if no == 3:
                no_data = int(line.split(",")[0].split()[1])
                dt = float(line.split(",")[1].split()[1])
            if no > 3:
                [acc.append(float(item)) for item in line.split()]
    time = [dt * item for item in range(no_data)]
    return(time, acc, dt)


def butter_bandpass(lowcut, highcut, fs, order=4):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a
def butter_bandpass_filter(data, lowcut, highcut, fs, order=2):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)

    y = signal.detrend( y , type='linear')

    return y


def target_spectrum(Ss, S1, soil):
    
    Ss_range = [0.25 , 0.50 , 0.75, 1.00 , 1.25 , 1.50 ]
    FS_table = {"ZA": [0.8 , 0.8 , 0.8 , 0.8 , 0.8 , 0.8], 
                "ZB": [0.9 , 0.9 , 0.9 , 0.9 , 0.9 , 0.9], 
                "ZC": [1.3 , 1.3 , 1.2 , 1.2 , 1.2 , 1.2],
                "ZD": [1.6 , 1.4 , 1.2 , 1.1 , 1.0 , 1.0],
                "ZE": [2.4 , 1.7 , 1.3 , 1.1 , 0.9 , 0.8]}

    S1_range = [0.10 , 0.20 , 0.30, 0.40 , 0.50 , 0.60 ]
    F1_table = {"ZA": [0.8 , 0.8 , 0.8 , 0.8 , 0.8 , 0.8], 
                "ZB": [0.8 , 0.8 , 0.8 , 0.8 , 0.8 , 0.8], 
                "ZC": [1.5 , 1.5 , 1.5 , 1.5 , 1.5 , 1.4],
                "ZD": [2.4 , 2.2 , 2.0 , 1.9 , 1.8 , 1.7],
                "ZE": [4.2 , 3.3 , 2.8 , 2.4 , 2.2 , 2.0]}


    if Ss < Ss_range[0]:
        FS_satir = np.polyfit(Ss_range[0:2], list(FS_table[soil])[0:2], 1)
        FS_katsayisi = np.poly1d( FS_satir )
        Fs = float( format(FS_katsayisi(Ss) , '.2f') )
        SDs = Ss * Fs
    elif Ss > Ss_range[-1]:
        FS_satir = np.polyfit(Ss_range[-3:-1], list(FS_table[soil])[-3:-1], 1)
        FS_katsayisi = np.poly1d( FS_satir )
        Fs = float( format(FS_katsayisi(Ss) , '.2f') )
        SDs = Ss * Fs    
    else:
        FS_satir = interp1d(Ss_range, FS_table[soil], kind='linear')
        FS_katsayisi = FS_satir(Ss)
        Fs = round( float(FS_katsayisi) , 2) 
        SDs = Ss * Fs

        
    if S1 < S1_range[0] :
        F1_satir = np.polyfit(S1_range[0:2], list(F1_table[soil])[0:2], 1)
        F1_katsayisi = np.poly1d( F1_satir )
        F1 = float( format(F1_katsayisi(S1) , '.2f') )
        SD1 = S1 * F1
    elif S1 > S1_range[-1]:
        F1_satir = np.polyfit(S1_range[-3:-1], list(F1_table[soil])[-3:-1], 1)
        F1_katsayisi = np.poly1d( F1_satir )
        F1 = float( format(F1_katsayisi(S1) , '.2f') )
        SD1 = S1 * F1

    else:    
        F1_satir = interp1d(S1_range, F1_table[soil], kind='linear')
        F1_katsayisi = F1_satir(S1)
        F1 = round(float(F1_katsayisi) , 2)
        SD1 = S1 * F1
        
    TA = 0.2 * SD1 / SDs
    TB = SD1 / SDs
    TL = 6
    
    T_range = np.arange(0, 6.05, 0.076).tolist()
    T_list = []
    for i in T_range:
        T_list.append(round(i,3))
        
    Sa = []
    
    for i in T_list:
        
        if i <TA:
            Sa.append(round((0.4 + 0.6*(i/TA))*SDs, 4))
            
        elif i >= TA and i<=TB:
            Sa.append(round(SDs, 4))
            
        elif i>TB and i <=TL:
            Sa.append(round(SD1/i, 4))
            
        elif i>TL:
            Sa.append(round(SD1*TL/(i**2), 4))
            
    target_spec = {"T" : T_list,
                   "Sa" : Sa}

    target_spec_df = pd.DataFrame().from_dict(target_spec)
    
    return target_spec_df


def SingleComponentScale(file_list, target_spectra, T1):

    ## Creating Response Spectras of the Records
    spectra_df = pd.DataFrame()

    for file_name in file_list:

        time, acceleration, dt = ReadPeer(os.path.join("./0-Data" , file_name))

        lowcut = .05
        highcut = 20.0

        acceleration_filtered = butter_bandpass_filter(acceleration, lowcut, highcut, 1 / dt, order=2) 
        acceleration_filtered = signal.detrend(acceleration_filtered, type='linear')

        T, Sd, Sv, Sa = ResponseSpectra(acceleration_filtered , sampling_interval = dt , damping_ratio = 0.05 , Tmax = 4, T_int= 0.05)

        spectra_df["Period"] = T
        spectra_df[ file_name.replace(".at2","")] = Sa


    plt.figure()
    carpim = {}
    for column in spectra_df.columns[1:] :
        plt.plot( spectra_df["Period"] , spectra_df[ column] ,  "gray", lw= 0.5 )
        carpim[ column ] = round( sum( ( target_spectra["Sa"]  * spectra_df[ column]) ) / sum( spectra_df[ column]**2 ) , 3)

    plt.plot( target_spectra["T"] , target_spectra["Sa"])

    plt.figure()

    factored_sa = pd.DataFrame()
    for key in carpim.keys() : 
        factored_sa[ key ] = list( carpim[ key ] * spectra_df[ key] )
        plt.plot( spectra_df["Period"] ,  factored_sa[ key ] , "gray", lw= 0.1)


    ortalama = factored_sa.mean( axis = 1 )

    T1_alt , T1_ust = round( 0.2 * T1  ,2 ) , round(1.5*T1  , 2)

    [ plt.axvline( item  , color= "k", linestyle = "-." ) for item in [T1 , T1_alt , T1_ust] ]

    indexAralik = spectra_df.index[(spectra_df["Period"] >= T1_alt) & (spectra_df["Period"] <= T1_ust)].tolist()

    plt.axvspan( T1_alt , T1_ust  , alpha=0.1, color='blue')


    plt.plot( spectra_df["Period"] , spectra_df.iloc[:,1:],"c" , lw = 0.1 )

    plt.plot( spectra_df["Period"] , spectra_df.iloc[:,1:].mean(axis=1),"r--" , lw = 1 , label="Non-Scaled Average" )

    plt.plot( target_spectra["T"] , target_spectra["Sa"], "blue" , lw= "2" , label="Target Spectrum" )

    plt.plot( target_spectra["T"] , ortalama , "r-" ,lw= "1" , label="Preliminary Scale Average" )

    SF_ortalama , tol , counter  = 1 , 1 , 0

    while tol > 0.01:
    
        ortalama_Scaled = ortalama.iloc[indexAralik] * SF_ortalama 

        farklar = target_spectra["Sa"].iloc[ indexAralik ] - ortalama_Scaled 
    
        if max(farklar) > 0 : 
            SF_ortalama  += 0.01
    
        if max(farklar) < 0 :
            SF_ortalama  -= 0.01
    
        if max( farklar ) > 0 and  max(farklar) < 0.01 : 
            tol = 0.001 

        counter += 1
        if counter == 50: 
            tol = 0.001
        #print(f"{counter} and { max(farklar) } and {SF_ortalama}") 

    plt.plot( target_spectra["T"] , ortalama * SF_ortalama , "green" ,lw= "2" , label="Last Scaled Average" )

    #plt.grid( color='grey', linestyle='-', linewidth=.1 , axis = "both")
    plt.xlabel("Period (s)"),plt.ylabel("Sa (g)")
    plt.legend()
    plt.xlim(left= 0 ), plt.ylim( bottom = 0 )
    plt.title( f"Scaled Spectra for T={T1}" )
    plt.tight_layout()
    plt.show()

    ScaleFactors = {}
    for key in carpim.keys():
        ScaleFactors[ key ] = round( carpim[key] * SF_ortalama , 3) 

    print( ScaleFactors )