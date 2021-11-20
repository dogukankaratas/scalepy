from numpy.lib.function_base import average
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pandas._config.config import reset_option
from pandas.core.frame import DataFrame
from scipy.interpolate import interp1d


def targetSpectrum(Ss, S1, soil):
    
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
    
    x_spectra = pd.read_csv("spectral_x.csv")
    cols = x_spectra.columns.tolist()
    t = []

    for i in cols[2:]:
        t.append(float(i))
    T_list = t
        
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


def recordSelection(magnitude_range: str = '4 9',
                   vs_range: str = '0 250', 
                   rjb_range: str = '0 250',
                   fault_mechnanism: str = 'Strike - Slip', 
                   duration_5_75_range: str = '0 30', 
                   duration_5_95_range: str = '0 50', 
                   arias_intensity_range: str = '0 5', 
                   target_spectrum: DataFrame = pd.DataFrame()):
    """
    Args:
        magnitude_range (str): Magnitude Range
        vs_range (str): VS30 Range
        rjb_range (str): RJB(km) Range
        fault_mechnanism (str): Fault Mechanism
            - Normal
            - Strike-Slip
            - Reverse
            - Reverse Oblique
            - Normal Oblique
        duration_5_75_range (str): 5-75% Duration(sec) Range
        duration_5_95_range (str): 5-95% Duration(sec) Range
        arias_intensity_range (str): Arias Intensity (m/sec) Range
        target_spectra (dataframe): Target Spectra Dataframe
    """

    # Read the Meta Data
    eqe_df = pd.read_csv("meta_data.csv")

    # Split Inputs
    min_m, max_m = [float(x) for x in magnitude_range.split()]
    min_vs, max_vs = [float(x) for x in vs_range.split()]
    min_r, max_r = [float(x) for x in rjb_range.split()]
    min_d_75, max_d_75 = [float(x) for x in duration_5_75_range.split()]
    min_d_95, max_d_95 = [float(x) for x in duration_5_95_range.split()]
    min_arias, max_arias = [float(x) for x in arias_intensity_range.split()]

    # Filter the Dataframe acc. to the Inputs
    eqe_s = eqe_df[(eqe_df[" Magnitude"] > min_m) & (eqe_df[" Magnitude"] < max_m) 
                      & (eqe_df[" Vs30 (m/sec)"] > min_vs) & (eqe_df[" Vs30 (m/sec)"] < max_vs)
                      &  (eqe_df[" Rjb (km)"] > min_r) & (eqe_df[" Rjb (km)"] < max_r)
                      &  (eqe_df[" 5-75% Duration (sec)"] > min_d_75) & (eqe_df[" 5-75% Duration (sec)"] < max_d_75)
                      &  (eqe_df[" 5-95% Duration (sec)"] > min_d_95) & (eqe_df[" 5-95% Duration (sec)"] < max_d_95)
                      &  (eqe_df[" Arias Intensity (m/sec)"] > min_arias) & (eqe_df[" Arias Intensity (m/sec)"] < max_arias)]


    if fault_mechnanism == 'Strike - Slip':
        eqe_s = eqe_s[(eqe_s[' Mechanism'] == ' strike slip')]
    
    elif fault_mechnanism ==  'Normal' :
        eqe_s = eqe_s[(eqe_s[' Mechanism'] == ' Normal')]
          
    elif fault_mechnanism ==  'Reverse' :
        eqe_s = eqe_s[(eqe_s[' Mechanism'] == ' Reverse')]

    elif fault_mechnanism ==  'Reverse - Oblique' :
        eqe_s = eqe_s[(eqe_s[' Mechanism'] == ' Reverse Oblique')]

    elif fault_mechnanism ==  'Normal - Oblique' :
        eqe_s = eqe_s[(eqe_s[' Mechanism'] == ' Normal Oblique')]

    else : 
        print("Invalid Mechanism!")

    x_spectra = pd.read_csv("spectral_x.csv")
    y_spectra = pd.read_csv("spectral_y.csv")

    rsn_selected = eqe_s[' Record Sequence Number'].tolist()

    x_df = pd.DataFrame()
    y_df = pd.DataFrame()

    geo_mean_df = pd.DataFrame()

    for i in rsn_selected:
        x_df[str(i)] = x_spectra.loc[x_spectra['RSN'] == i].iloc[0].tolist()[2:]
        y_df[str(i)] = y_spectra.loc[y_spectra['RSN'] == i].iloc[0].tolist()[2:]

    # Select 11 records which has miminum difference between target spectra
    for i in rsn_selected:
        geo_mean_df[str(i)] = [(x*y)**(1/2) for x,y in zip(x_df[str(i)].to_list(), y_df[str(i)].to_list())]

    differ_dict = {}

    sa_list = target_spectrum['Sa'].tolist()
    for i in geo_mean_df.columns:
        differ_dict[i] = [abs(x-y) for x,y in zip(sa_list, geo_mean_df[str(i)])]
    for i in geo_mean_df.columns:
        differ_dict[i] = max(differ_dict[i])
    key_list = []
    for key, value in differ_dict.items():
        if value in sorted(differ_dict.values())[:11]:
            key_list.append(key)

    return key_list

def amplitudeScaling(keys, target, period):

    spectral_data_x = pd.read_csv("spectral_x.csv")
    spectral_data_y = pd.read_csv("spectral_y.csv")

    key_int = []
    for i in keys:
        key_int.append(int(i))
    selected_x = spectral_data_x.loc[spectral_data_x['RSN'].isin(key_int)]
    selected_y = spectral_data_y.loc[spectral_data_y['RSN'].isin(key_int)]

    t_str = selected_x.columns.tolist()[2:]
    t = []
    for i in t_str:
        t.append(float(i))
    rsn_selected = selected_x['RSN'].tolist()

    eqe_selected_x = pd.DataFrame()
    eqe_selected_x.insert(0, 'T', t)
    eqe_selected_y = pd.DataFrame()
    
    eqe_selected_y.insert(0, 'T', t)

    for i in rsn_selected:
        eqe_selected_x[str(i)] = selected_x.loc[selected_x['RSN'] == i].iloc[0].tolist()[2:]
        eqe_selected_y[str(i)] = selected_y.loc[selected_y['RSN'] == i].iloc[0].tolist()[2:]

    def geomean_func(acc_1 , acc_2):
        geo_mean = []
        for i , j in zip(acc_1 , acc_2):
            geo_mean.append( round( ( i * j)**(0.5) , 4 ) )
        return( geo_mean)

    def srss_func(acc_1 , acc_2):
        srss_mean = []
        for i , j in zip(acc_1 , acc_2):
            srss_mean.append( round( ( i**2 + j**2)**(0.5) , 4 ) )
        return( srss_mean)

    geo_mean_df = pd.DataFrame()
    geo_mean_df.insert(0, 'T', t)

    for i in rsn_selected:
        geo_mean_df[str(i)] = geomean_func(eqe_selected_x[str(i)].tolist(), eqe_selected_y[str(i)].tolist())
        
    filtered_target = target[(target["T"] > 0.2*int(period)) & (target["T"] < 1.5*int(period))]
    filtered_geo_mean = geo_mean_df[(geo_mean_df["T"] > 0.2*int(period)) & (geo_mean_df["T"] < 1.5*int(period))]

    geo_sf_dict = {}
    for i in filtered_geo_mean.columns[1:]:
        num = 0
        denom = 0
        for x,y in zip(filtered_geo_mean[str(i)].tolist(), filtered_target['Sa'].tolist()):
            num += x*y
            denom += y**2
            geo_sf = (num/denom)
        geo_sf_dict[i] = geo_sf

    multipilied_selected_x = eqe_selected_x.copy()
    multipilied_selected_y = eqe_selected_y.copy()

    for i in rsn_selected:
        multipilied_selected_x[str(i)] = geo_sf_dict[str(i)] * multipilied_selected_x[str(i)]
        multipilied_selected_y[str(i)] = geo_sf_dict[str(i)] * multipilied_selected_y[str(i)]

    srss_mean_df = pd.DataFrame()
    srss_mean_df.insert(0, 'T', t)

    for i in rsn_selected:
        srss_mean_df[str(i)] = srss_func(multipilied_selected_x[str(i)].tolist(), multipilied_selected_y[str(i)].tolist())

    rsn_str = []
    for i in rsn_selected:
        rsn_str.append(str(i))
    srss_mean_df['Mean'] = srss_mean_df[rsn_str].mean(axis=1)

    filtered_srss= srss_mean_df[(srss_mean_df["T"] > 0.2*int(period)) & (srss_mean_df["T"] < 1.5*int(period))]

    num = 0
    denom = 0
    for x,y in zip(filtered_srss['Mean'].tolist(), filtered_target['Sa'].tolist()):
        num += x*y
        denom += y**2
    srss_sf = (num/denom)
    
    final_sf_dict = {}
    for i,j in zip(geo_sf_dict.keys(), geo_sf_dict.values()):
        final_sf_dict[i] = srss_sf * j
    
    print(final_sf_dict)

    srss_mean_df_final = srss_mean_df.copy()

    for i,j in zip(rsn_selected, final_sf_dict.values()):
        srss_mean_df_final[str(i)] = srss_mean_df_final[str(i)] * j

    srss_mean_df_final['Mean'] = srss_mean_df_final[rsn_str].mean(axis=1)
    
    print(srss_mean_df_final)

    


    # Visualization
    for i,j in zip(rsn_selected, final_sf_dict.values()):
        plt.plot(eqe_selected_x['T'], j * eqe_selected_x[str(i)], linewidth =0.5, color = 'grey')
        plt.plot(eqe_selected_y['T'], j * eqe_selected_y[str(i)], linewidth =0.5, color = 'grey')



    plt.plot(srss_mean_df_final['T'], srss_mean_df_final['Mean'], color = 'blue')
    plt.xlim(0, 3)
    plt.plot(target['T'], target['Sa'], color = 'red', linewidth = 3)


