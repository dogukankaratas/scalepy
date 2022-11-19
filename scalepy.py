import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pandas.core.frame import DataFrame
from scipy.interpolate import interp1d
import similaritymeasures
from statistics import median

# Ignores performance warnings
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
simplefilter(action="ignore", category=FutureWarning)

def targetSpectrum(Ss, S1, soil):
    """
    Args:
       Ss: Spectral Acceleration Parameter at Short Periods
       S1: Spectral Acceleration Parameter at 1-sec
       soil: Soil Type
    """
    
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
    
    x_spectra = pd.read_csv("data/spectral_x.csv")
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
                   target_spectrum: DataFrame = pd.DataFrame() , 
                   pulse_type : str = "Any",
                   period: float = 1):
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
        pulse (int): Pulse [Pulse‐like (1)] or Non-pulse [non‐pulse‐like (0)] or Any[any (2)] indicator
    """
    print( "-" * 50 )
    for key , val in {"magnitude_range":magnitude_range , "vs_range" : vs_range, "rjb_range" :rjb_range, "fault_mechnanism": fault_mechnanism, "duration_5_75_range" : duration_5_75_range,"duration_5_95_range" : duration_5_95_range, "arias_intensity_range" : arias_intensity_range, "target_spectrum" : 1, "pulse_type" : pulse_type ,"period" :period}.items() :
        print( f"{key.casefold()} = {val}")

    print( "-" * 50 )

    # Read the Meta Data
    eqe_df = pd.read_csv("data/meta_data-R1.csv")

    print( f"Total number of PEER EQE Record is = {len( eqe_df)}")

    # Split Inputs
    min_m, max_m = [float(x) for x in magnitude_range.split()]
    min_vs, max_vs = [float(x) for x in vs_range.split()]
    min_r, max_r = [float(x) for x in rjb_range.split()]
    min_d_75, max_d_75 = [float(x) for x in duration_5_75_range.split()]
    min_d_95, max_d_95 = [float(x) for x in duration_5_95_range.split()]
    min_arias, max_arias = [float(x) for x in arias_intensity_range.split()]

    # Filter the Dataframe acc. to the Inputs
    eqe_s = eqe_df[(eqe_df["Magnitude"] >= min_m) & (eqe_df["Magnitude"] <= max_m) 
                      & (eqe_df["Vs30(m/sec)"] >= min_vs) & (eqe_df["Vs30(m/sec)"] <= max_vs)
                      &  (eqe_df["Rjb(km)"] >= min_r) & (eqe_df["Rjb(km)"] <= max_r)
                      &  (eqe_df["5-75%Duration(sec)"] >= min_d_75) & (eqe_df["5-75%Duration(sec)"] <= max_d_75)
                      &  (eqe_df["5-95%Duration(sec)"] >= min_d_95) & (eqe_df["5-95%Duration(sec)"] <= max_d_95)
                      &  (eqe_df["AriasIntensity(m/sec)"] >= min_arias) & (eqe_df["AriasIntensity(m/sec)"] <= max_arias)]
    
    # Pulse type filtering
    if pulse_type == "Pulse" :
        eqe_s = eqe_s[ eqe_s["Pulse"] == 1 ]
    elif pulse_type == "Non-Pulse" :
        eqe_s = eqe_s[ eqe_s["Pulse"] == 2 ]
    else :
        pass

    # Mechanism type filtering
    if fault_mechnanism == 'Strike - Slip':
        eqe_s = eqe_s[(eqe_s['Mechanism'] == ' strike slip')]
    
    elif fault_mechnanism ==  'Normal' :
        eqe_s = eqe_s[(eqe_s['Mechanism'] == ' Normal')]
          
    elif fault_mechnanism ==  'Reverse' :
        eqe_s = eqe_s[(eqe_s['Mechanism'] == ' Reverse')]

    elif fault_mechnanism ==  'Reverse - Oblique' :
        eqe_s = eqe_s[(eqe_s['Mechanism'] == ' Reverse Oblique')]

    elif fault_mechnanism ==  'Normal - Oblique' :
        eqe_s = eqe_s[(eqe_s['Mechanism'] == ' Normal Oblique')]

    elif fault_mechnanism == 'Oblique':
        eqe_s = eqe_s[((eqe_s['Mechanism'] == ' Normal Oblique') | (eqe_s['Mechanism'] == ' Reverse Oblique'))]

    else : 
        print("Invalid Mechanism!")
    
    print( f"Total number of appropriate PEER EQE Record (1st selection) is = {len( eqe_s)}")

    # Select 11 records with minimal difference and 3 from same earthquake
    first_three_list = []

    for counter, temp_df in eqe_s.groupby( by = "EarthquakeName") :
        [first_three_list.append(item) for item in  temp_df["RecordSequenceNumber"][:3].to_list() ]
        
    print( f"Total number of appropriate PEER EQE Record (2nd selection) is = {len( first_three_list)}")

    eqe_s_filtered = eqe_s[ eqe_s["RecordSequenceNumber"].isin( first_three_list )]

    # Selected Record Sequence Numbers
    rsn_selected = eqe_s_filtered['RecordSequenceNumber'].tolist()

    # Read Spectral Data
    spectral_data_x = pd.read_csv("data/spectral_x.csv")
    spectral_data_y = pd.read_csv("data/spectral_y.csv")

    spectra_selected_x = spectral_data_x.loc[spectral_data_x['RSN'].isin(rsn_selected)]
    spectra_selected_y = spectral_data_y.loc[spectral_data_y['RSN'].isin(rsn_selected)]

    t_str = spectra_selected_x.columns.tolist()[2:]
    t = []
    for i in t_str:
        t.append(float(i))

    # Reverse the DataFrames
    eqe_selected_x = pd.DataFrame()
    eqe_selected_x.insert(0, 'T', t)
    eqe_selected_y = pd.DataFrame()
    eqe_selected_y.insert(0, 'T', t)

    for i in rsn_selected:
        eqe_selected_x[ i ] = spectra_selected_x.loc[spectra_selected_x['RSN'] == i].iloc[0].tolist()[2:]
        eqe_selected_y[ i ] = spectra_selected_y.loc[spectra_selected_y['RSN'] == i].iloc[0].tolist()[2:]

    # Calculation of Geometric Mean of the Records
    geo_mean_df = pd.DataFrame()
    geo_mean_df['T'] = t

    for i in rsn_selected:
        geo_mean_df[i ] = [(x*y)**(1/2) for x,y in zip(eqe_selected_x[ i ].to_list(), eqe_selected_y[ i ].to_list())]

    # range of interest
    geo_mean_df_range = geo_mean_df[ ( geo_mean_df["T"] >= 0.2 * period) & ( geo_mean_df["T"] <= 1.5 * period ) ]
    target_range = target_spectrum[ ( target_spectrum["T"] >= 0.2 * period) & ( target_spectrum["T"] <= 1.5 * period )]

    # Create 2D Arrays for Records
    similarities_df = pd.DataFrame()
    similarities_df['RSN'] = rsn_selected
    df = []
    area = []
    cl = []
    dtw = []
    mae = []
    mse = []

    target_array = np.column_stack((target_range['T'], target_range['Sa']))

    for i in rsn_selected:
        df.append(similaritymeasures.frechet_dist(np.column_stack((geo_mean_df_range['T'], geo_mean_df_range[i])), target_array))
        area.append(similaritymeasures.area_between_two_curves(np.column_stack((geo_mean_df_range['T'], geo_mean_df_range[i])), target_array))
        cl.append(similaritymeasures.curve_length_measure(np.column_stack((geo_mean_df_range['T'], geo_mean_df_range[i])), target_array))
        dtw.append(similaritymeasures.dtw(np.column_stack((geo_mean_df_range['T'], geo_mean_df_range[i])), target_array)[0])
        mae.append(similaritymeasures.mae(np.column_stack((geo_mean_df_range['T'], geo_mean_df_range[i])), target_array))
        mse.append(similaritymeasures.mse(np.column_stack((geo_mean_df_range['T'], geo_mean_df_range[i])), target_array))

    similarities_df['DF'] = df
    similarities_df['AREA'] = area
    similarities_df['CL'] = cl
    similarities_df['DTW'] = dtw
    similarities_df['MAE'] = mae
    similarities_df['MSE'] = mse

    similarityMean = []
    for d,a,c,w,e,s in zip(df, area, cl, dtw, mae, mse):
        similarityMean.append((d+a+c+w+e+s)/6)

    similarities_df['Mean'] = similarityMean

    # get most similiar 11 records
    similarities_df = similarities_df.sort_values(by='Mean')
    selected_keys = similarities_df.head(11)['RSN'].to_list()

    # all the available records
    plt.figure(figsize= (10, 5))
    for eqe_name in rsn_selected : 
        plt.plot( t , eqe_selected_x[eqe_name] , linestyle = "--" , lw = .7 , color = "gray")
        plt.plot( t , eqe_selected_y[eqe_name] , linestyle = "-." , lw = .7 , color = "gray")
    plt.plot( t , target_spectrum["Sa"])
    plt.axvspan( 0.2 * float(period) , 1.5 * float( period) , facecolor = "green" , alpha = 0.1)
    plt.xlim( 0  , 3 )
    plt.ylim(bottom = 0 )
    plt.box(False)
    plt.axhline(c="k")
    plt.axvline(c="k")
    plt.title("Unoptimized Filtered Records")

    # optimum selected records visualization
    plt.figure(figsize=(10,5))
    for i in selected_keys:
        plt.plot(t, eqe_selected_x[i], linestyle = "--" , lw = .7 , color = "gray")
        plt.plot(t, eqe_selected_y[i], linestyle = "--" , lw = .7 , color = "gray")
    plt.plot(t, target_spectrum['Sa'])
    plt.axvspan( 0.2 * float( period) , 1.5 * float( period) , facecolor = "green" , alpha = 0.1)
    plt.xlim(0, 3)
    plt.ylim(bottom = 0 )
    plt.box(False)
    plt.axhline(c="k")
    plt.axvline(c="k")
    plt.title("Optimum Selected Records")

    return selected_keys


    ## Function for Amplitude Scaling ##

def amplitudeScaling(key_list, target , period, components = 'srss'):

    eqe_df = pd.read_csv("data/meta_data-R1.csv")

    spectral_data_x = pd.read_csv("data/spectral_x.csv")
    spectral_data_y = pd.read_csv("data/spectral_y.csv")

    selected_x = spectral_data_x.loc[spectral_data_x['RSN'].isin(key_list)]
    selected_y = spectral_data_y.loc[spectral_data_y['RSN'].isin(key_list)]

    t_str = selected_x.columns.tolist()[2:]
    t = []
    for i in t_str:
        t.append(float(i))

    rsn_selected = key_list

    eqe_selected_x = pd.DataFrame()
    eqe_selected_x.insert(0, 'T', t)
    eqe_selected_y = pd.DataFrame()
    eqe_selected_y.insert(0, 'T', t)

    for i in rsn_selected:
        eqe_selected_x[ i ] = selected_x.loc[selected_x['RSN'] == i].iloc[0].tolist()[2:]
        eqe_selected_y[ i ] = selected_y.loc[selected_y['RSN'] == i].iloc[0].tolist()[2:]

    # Visualization Before Scaling
    plt.figure(figsize= (10, 5) )
    for eqe_name in key_list : 
        plt.plot( t , eqe_selected_x[eqe_name] , linestyle = "--" , lw = .7 , color = "gray")
        plt.plot( t , eqe_selected_y[eqe_name] , linestyle = "-." , lw = .7 , color = "gray")
    plt.plot( t , target["Sa"])
    plt.axvspan( 0.2 * float( period) , 1.5 * float( period) , facecolor = "green" , alpha = 0.1)
    plt.xlim( 0  , 3 )
    plt.ylim(bottom = 0 )
    plt.box(False)
    plt.axhline(c="k")
    plt.axvline(c="k")
    plt.title("Unscaled Spectra of the selected ground motions")

    # Create Geometric Mean, SRSS Mean, RotD50 and RotD100 Functions
    def geomean_func(acc_1 , acc_2):
        geo_mean = []
        for i , j in zip(acc_1 , acc_2):
            geo_mean.append( round( ( i * j)**(0.5) , 4 ) )

        return geo_mean

    def srss_func(acc_1 , acc_2):
        srss_mean = []
        for i , j in zip(acc_1 , acc_2):
            srss_mean.append( round( ( i**2 + j**2)**(0.5) , 4 ) )

        return srss_mean

    def rotD50_func(acc_1, acc_2):
        rot50 = []
        for i,j in zip(acc_1, acc_2):
            rot50.append(np.percentile(np.array([i, j]), 50))

        return rot50

    def rotD100_func(acc_1, acc_2):
        rot100 = []
        for i,j in zip(acc_1, acc_2):
            rot100.append(np.percentile(np.array([i, j]), 100))
        
        return rot100

    # Find Geometric Mean
    geo_mean_df = pd.DataFrame()
    geo_mean_df.insert(0, 'T', t)

    for i in rsn_selected:
        geo_mean_df[  i ] = geomean_func( eqe_selected_x[ i ].tolist() , eqe_selected_y[ i ].tolist() )
        
    plt.figure( figsize = (10 ,5))
    color_list = ("blue", "orange", "green", "red", "purple", "brown", "pink", "gray", "olive", "cyan", "lime")
    for rsn, color in zip(key_list, color_list) : 
        plt.plot( t , geo_mean_df[ rsn] , c = color , label = rsn)
    plt.plot( t , target["Sa"] , c = "red")
    plt.legend( loc = 1)
    plt.axvspan( 0.2 * float( period) , 1.5 * float( period) , facecolor = "green" , alpha = 0.1)
    plt.xlim( 0  , 3 )
    plt.ylim(bottom = 0 )
    plt.box(False)
    plt.axhline(c="k")
    plt.axvline(c="k")
    plt.title("Geomean Spectra of the selected ground motions")       
    
    # Slice the period range from target and geomean spectra
    filtered_target = target[(target["T"] >= 0.2*float(period)) & (target["T"] <= 1.5*float(period))]

    filtered_geo_mean = geo_mean_df[(geo_mean_df["T"] >= 0.2*float(period)) & (geo_mean_df["T"] <= 1.5*float(period))]

    # Find the differences in period range 
    geo_sf_dict = {}
    num = 0
    denom = 0
    for rsn in key_list:
        for x,y in zip(filtered_geo_mean[ rsn ].to_list(), filtered_target[ "Sa" ].to_list()):
            num += x*y
            denom += x**2
        geo_sf_dict[rsn] = (num/denom)

    # First scaling= Geomen x geo_sf
    multiplied_selected_x = eqe_selected_x.copy()
    multiplied_selected_y = eqe_selected_y.copy()

    for rsn in rsn_selected:
        multiplied_selected_x[ rsn ] = geo_sf_dict[ rsn ] * multiplied_selected_x[ rsn ]
        multiplied_selected_y[ rsn ] = geo_sf_dict[ rsn ] * multiplied_selected_y[ rsn ]

    geo_mean_1st_scaled_df = pd.DataFrame()
    geo_mean_1st_scaled_df["T"] = geo_mean_df["T"]
    for rsn in rsn_selected : 
        geo_mean_1st_scaled_df[ rsn ] = geo_sf_dict[ rsn] *  geo_mean_df[ rsn]
    geo_mean_1st_scaled_df[ "Mean" ]  = geo_mean_1st_scaled_df[rsn_selected].mean( axis = 1 )

    plt.figure( figsize = (10 ,5))
    for rsn in key_list : 
        plt.plot( t , multiplied_selected_x[ rsn] , c = "gray" , lw = 0.5 , label = rsn)
        plt.plot( t , multiplied_selected_y[ rsn] , c = "gray" , lw = 0.5 , label = rsn)
    plt.plot( t , target["Sa"] , c = "red")
    plt.plot( t , geo_mean_1st_scaled_df[ "Mean"] , c = "blue" , label = "Geomean-Scaled-Average")

    plt.axvspan( 0.2 * float( period) , 1.5 * float( period) , facecolor = "green" , alpha = 0.1)
    plt.xlim( 0  , 3 )
    plt.ylim(bottom = 0 )
    plt.box(False)
    plt.axhline(c="k")
    plt.axvline(c="k")
    plt.title("Average Geomean Spectra of the scaled selected ground motions")     

    #################################### Use SRSS Function To Find Spectral Component Mean ####################################  

    if components == 'srss':    
        srss_mean_df = pd.DataFrame()
        srss_mean_df.insert(0, 'T', t)
        for rsn in rsn_selected:
            srss_mean_df[  rsn ] = srss_func(multiplied_selected_x[ rsn ].tolist(), multiplied_selected_y[ rsn ].tolist())
        rsn_str = []
        for i in rsn_selected:
            rsn_str.append( str(i))
                    
        srss_mean_df['Mean'] = srss_mean_df[ key_list].mean(axis=1)
        plt.figure( figsize = (10 ,5))
        for rsn in key_list : 
            plt.plot( t , multiplied_selected_x[ rsn] , c = "gray" , lw = 0.5 , label = rsn)
            plt.plot( t , multiplied_selected_y[ rsn] , c = "gray" , lw = 0.5 , label = rsn)
        plt.plot( t , target["Sa"] , c = "red")
        plt.plot( t , geo_mean_1st_scaled_df[ "Mean"] , c = "blue" , label = "Geomean-Scaled")
        plt.plot( t , srss_mean_df[ "Mean"] , c = "black" , label = "Srss-Scaled")
        plt.axvspan( 0.2 * float( period) , 1.5 * float( period) , facecolor = "green" , alpha = 0.1)
        plt.xlim( 0  , 3 )
        plt.ylim(bottom = 0 )
        plt.box(False)
        plt.axhline(c="k")
        plt.axvline(c="k")
        plt.title("Geomean and SRSS Spectra of the 1st scaled selected ground motions")     
        
        filtered_srss = srss_mean_df[(srss_mean_df["T"] >= 0.2*float(period)) & (srss_mean_df["T"] <= 1.5*float(period))]
        
        SF_ortalama , tol , counter  = 1 , 1 , 0
        spectra_increament_factor = 1.3
        while tol > 0.01:
            
            ortalama_Scaled = filtered_srss["Mean"] * SF_ortalama 
            farklar = spectra_increament_factor * filtered_target["Sa"] - ortalama_Scaled 
            
            if max(farklar) > 0 : 
                SF_ortalama  += 0.01
            
            if max(farklar) < 0 :
                SF_ortalama  -= 0.01
            
            if max( farklar ) > 0 and  max(farklar) < 0.01 : 
                tol = 0.001 
            counter += 1
            if counter == 50: 
                tol = 0.001

        #Spectral Increament Double Check
        increased_target = spectra_increament_factor * filtered_target["Sa"]
        differences = {}
        index = 0
        for i,j in zip(increased_target.to_list(), ortalama_Scaled.to_list()):
            differences[index] = round((i/j), 4)        
            index +=1
        max_dif = max(differences.values())
        index_max = [key for key, j in differences.items() if j == max_dif]
        inc = increased_target.to_list()[index_max[0]]/ortalama_Scaled.to_list()[index_max[0]]
        
        #Obtain Scale Factor Dictionary
        sf_dict = {}
        for key , val in geo_sf_dict.items(): 
            sf_dict[ key ] = round( SF_ortalama * val * inc , 4 )

        
        # Obtain the scaled spectral values
        eqe_selected_scaled_x , eqe_selected_scaled_y = pd.DataFrame() , pd.DataFrame()
        srss_mean_scaled_df = pd.DataFrame()
        srss_mean_scaled_df["T"] = t 
        for rsn in rsn_selected :
            eqe_selected_scaled_x[ rsn ] = sf_dict[ rsn ] * eqe_selected_x[ rsn]
            eqe_selected_scaled_y[ rsn ] = sf_dict[ rsn ] * eqe_selected_y[ rsn]
            srss_mean_scaled_df[  rsn ] = srss_func( eqe_selected_scaled_x[ rsn ].tolist(), eqe_selected_scaled_y[ rsn ].tolist())
        
        srss_mean_scaled_df["Mean"] = srss_mean_scaled_df[ rsn_selected ].mean( axis = 1 )
            
        plt.figure( figsize = (10 ,5))
        for rsn in key_list : 
            plt.plot( t , multiplied_selected_x[ rsn] , c = "gray" , lw = 0.5)
            plt.plot( t , multiplied_selected_y[ rsn] , c = "gray" , lw = 0.5)
        plt.plot( t , target["Sa"] , linestyle = "--" , c = "red", label= "Target")
        plt.plot( t , spectra_increament_factor *  target["Sa"] , c = "red" , lw = 1.5 ,  label = f"{spectra_increament_factor}x Target")
        plt.plot( t , geo_mean_1st_scaled_df[ "Mean"] , c = "blue" , label = "Geomean-Scaled")
        plt.plot( t , srss_mean_df[ "Mean"] , linestyle = "--" , c = "black" , label = "Srss-Scaled")
        plt.plot( t , srss_mean_scaled_df[ "Mean"] , c = "black" , lw = 1.5 , label = "SRSS-Mean-Scaled")
        plt.legend( loc = 1)
        plt.axvspan( 0.2 * float( period) , 1.5 * float( period) , facecolor = "green" , alpha = 0.1)
        plt.xlim( 0  , 3 )
        plt.ylim(bottom = 0 )
        plt.box(False)
        plt.axhline(c="k")
        plt.axvline(c="k")
        plt.title("Geomean and SRSS Spectra of the 2nd scaled selected ground motions")  ,
        plt.show()
        
        print( "-"*50)
        print( "Selected ground motions and scale factors")
        for rsn in rsn_selected : 
            temp_df = eqe_df[ eqe_df['RecordSequenceNumber'] == rsn]
            print( f"RSN{rsn} | { list( temp_df.EarthquakeName )[0] } | { sf_dict[ rsn ]}") 
        print( "-"*50)

            #################################### Use RotD50 Function To Find Spectral Component Mean ####################################

    elif components == 'rotd50':    
        rotd50_mean_df = pd.DataFrame()
        rotd50_mean_df.insert(0, 'T', t)
        for rsn in rsn_selected:
            rotd50_mean_df[  rsn ] = rotD50_func(multiplied_selected_x[ rsn ].tolist(), multiplied_selected_y[ rsn ].tolist())
        rsn_str = []
        for i in rsn_selected:
            rsn_str.append( str(i))
                    
        rotd50_mean_df['Mean'] = rotd50_mean_df[ key_list].mean(axis=1)
        plt.figure( figsize = (10 ,5))
        for rsn in key_list : 
            plt.plot( t , multiplied_selected_x[ rsn] , c = "gray" , lw = 0.5 , label = rsn)
            plt.plot( t , multiplied_selected_y[ rsn] , c = "gray" , lw = 0.5 , label = rsn)
        plt.plot( t , target["Sa"] , c = "red")
        plt.plot( t , geo_mean_1st_scaled_df[ "Mean"] , c = "blue" , label = "Geomean-Scaled")
        plt.plot( t , rotd50_mean_df[ "Mean"] , c = "black" , label = "RotD50-Scaled")
        plt.axvspan( 0.2 * float( period) , 1.5 * float( period) , facecolor = "green" , alpha = 0.1)
        plt.xlim( 0  , 3 )
        plt.ylim(bottom = 0 )
        plt.box(False)
        plt.axhline(c="k")
        plt.axvline(c="k")
        plt.title("Geomean and RotD50 Spectra of the 1st scaled selected ground motions")     
        
        filtered_rotd50 = rotd50_mean_df[(rotd50_mean_df["T"] >= 0.2*float(period)) & (rotd50_mean_df["T"] <= 1.5*float(period))]
        
        SF_ortalama , tol , counter  = 1 , 1 , 0
        spectra_increament_factor = 1.3
        while tol > 0.01:
            
            ortalama_Scaled = filtered_rotd50["Mean"] * SF_ortalama 
            farklar = spectra_increament_factor * filtered_target["Sa"] - ortalama_Scaled 
            
            if max(farklar) > 0 : 
                SF_ortalama  += 0.01
            
            if max(farklar) < 0 :
                SF_ortalama  -= 0.01
            
            if max( farklar ) > 0 and  max(farklar) < 0.01 : 
                tol = 0.001 
            counter += 1
            if counter == 50: 
                tol = 0.001

        #Spectral Increament Double Check
        increased_target = spectra_increament_factor * filtered_target["Sa"]
        differences = {}
        index = 0
        for i,j in zip(increased_target.to_list(), ortalama_Scaled.to_list()):
            differences[index] = round((i-j), 4)        
            index +=1
        max_dif = max(differences.values())
        index_max = [key for key, j in differences.items() if j == max_dif]
        inc = increased_target.to_list()[index_max[0]]/ortalama_Scaled.to_list()[index_max[0]]
        
        #Obtain Scale Factor Dictionary
        sf_dict = {}
        for key , val in geo_sf_dict.items(): 
            sf_dict[ key ] = round( SF_ortalama * val * inc , 4 )
        
        # Obtain the scaled spectral values
        eqe_selected_scaled_x , eqe_selected_scaled_y = pd.DataFrame() , pd.DataFrame()
        rotd50_mean_scaled_df = pd.DataFrame()
        rotd50_mean_scaled_df["T"] = t 
        for rsn in rsn_selected :
            eqe_selected_scaled_x[ rsn ] = sf_dict[ rsn ] * eqe_selected_x[ rsn]
            eqe_selected_scaled_y[ rsn ] = sf_dict[ rsn ] * eqe_selected_y[ rsn]
            rotd50_mean_scaled_df[  rsn ] = rotD50_func( eqe_selected_scaled_x[ rsn ].tolist(), eqe_selected_scaled_y[ rsn ].tolist())
        
        rotd50_mean_scaled_df["Mean"] = rotd50_mean_scaled_df[ rsn_selected ].mean( axis = 1 )
            
        plt.figure( figsize = (10 ,5))
        for rsn in key_list : 
            plt.plot( t , multiplied_selected_x[ rsn] , c = "gray" , lw = 0.5)
            plt.plot( t , multiplied_selected_y[ rsn] , c = "gray" , lw = 0.5)
        plt.plot( t , target["Sa"] , linestyle = "--" , c = "red", label= "Target")
        plt.plot( t , spectra_increament_factor *  target["Sa"] , c = "red" , lw = 1.5 ,  label = f"{spectra_increament_factor}x Target")
        plt.plot( t , geo_mean_1st_scaled_df[ "Mean"] , c = "blue" , label = "Geomean-Scaled")
        plt.plot( t , rotd50_mean_df[ "Mean"] , linestyle = "--" , c = "black" , label = "RotD50-Scaled")
        plt.plot( t , rotd50_mean_scaled_df[ "Mean"] , c = "black" , lw = 1.5 , label = "RotD50-Mean-Scaled")
        plt.legend( loc = 1)
        plt.axvspan( 0.2 * float( period) , 1.5 * float( period) , facecolor = "green" , alpha = 0.1)
        plt.xlim( 0  , 3 )
        plt.ylim(bottom = 0 )
        plt.box(False)
        plt.axhline(c="k")
        plt.axvline(c="k")
        plt.title("Geomean and RotD50 Spectra of the 2nd scaled selected ground motions")  ,
        plt.show()
        
        print( "-"*50)
        print( "Selected ground motions and scale factors")
        for rsn in rsn_selected : 
            temp_df = eqe_df[ eqe_df['RecordSequenceNumber'] == rsn]
            print( f"RSN{rsn} | { list( temp_df.EarthquakeName )[0] } | { sf_dict[ rsn ]}") 
        print( "-"*50)

    #################################### Use RotD100 Function To Find Spectral Component Mean ####################################

    elif components == 'rotd100':    
        rotd100_mean_df = pd.DataFrame()
        rotd100_mean_df.insert(0, 'T', t)
        for rsn in rsn_selected:
            rotd100_mean_df[  rsn ] = rotD100_func(multiplied_selected_x[ rsn ].tolist(), multiplied_selected_y[ rsn ].tolist())
        rsn_str = []
        for i in rsn_selected:
            rsn_str.append( str(i))
                    
        rotd100_mean_df['Mean'] = rotd100_mean_df[ key_list].mean(axis=1)
        plt.figure( figsize = (10 ,5))
        for rsn in key_list : 
            plt.plot( t , multiplied_selected_x[ rsn] , c = "gray" , lw = 0.5 , label = rsn)
            plt.plot( t , multiplied_selected_y[ rsn] , c = "gray" , lw = 0.5 , label = rsn)
        plt.plot( t , target["Sa"] , c = "red")
        plt.plot( t , geo_mean_1st_scaled_df[ "Mean"] , c = "blue" , label = "Geomean-Scaled")
        plt.plot( t , rotd100_mean_df[ "Mean"] , c = "black" , label = "RotD100-Scaled")
        plt.axvspan( 0.2 * float( period) , 1.5 * float( period) , facecolor = "green" , alpha = 0.1)
        plt.xlim( 0  , 3 )
        plt.ylim(bottom = 0 )
        plt.box(False)
        plt.axhline(c="k")
        plt.axvline(c="k")
        plt.title("Geomean and RotD100 Spectra of the 1st scaled selected ground motions")     
        
        filtered_rotd100 = rotd100_mean_df[(rotd100_mean_df["T"] >= 0.2*float(period)) & (rotd100_mean_df["T"] <= 1.5*float(period))]
        
        SF_ortalama , tol , counter  = 1 , 1 , 0
        spectra_increament_factor = 1.3
        while tol > 0.01:
            
            ortalama_Scaled = filtered_rotd100["Mean"] * SF_ortalama 
            farklar = spectra_increament_factor * filtered_target["Sa"] - ortalama_Scaled 
            
            if max(farklar) > 0 : 
                SF_ortalama  += 0.01
            
            if max(farklar) < 0 :
                SF_ortalama  -= 0.01
            
            if max( farklar ) > 0 and  max(farklar) < 0.01 : 
                tol = 0.001 
            counter += 1
            if counter == 50: 
                tol = 0.001

        #Spectral Increament Double Check
        increased_target = spectra_increament_factor * filtered_target["Sa"]
        differences = {}
        index = 0
        for i,j in zip(increased_target.to_list(), ortalama_Scaled.to_list()):
            differences[index] = round((i-j), 4)        
            index +=1
        max_dif = max(differences.values())
        index_max = [key for key, j in differences.items() if j == max_dif]
        inc = increased_target.to_list()[index_max[0]]/ortalama_Scaled.to_list()[index_max[0]]
        
        #Obtain Scale Factor Dictionary
        sf_dict = {}
        for key , val in geo_sf_dict.items(): 
            sf_dict[ key ] = round( SF_ortalama * val * inc , 4 )
        
        # Obtain the scaled spectral values
        eqe_selected_scaled_x , eqe_selected_scaled_y = pd.DataFrame() , pd.DataFrame()
        rotd100_mean_scaled_df = pd.DataFrame()
        rotd100_mean_scaled_df["T"] = t 
        for rsn in rsn_selected :
            eqe_selected_scaled_x[ rsn ] = sf_dict[ rsn ] * eqe_selected_x[ rsn]
            eqe_selected_scaled_y[ rsn ] = sf_dict[ rsn ] * eqe_selected_y[ rsn]
            rotd100_mean_scaled_df[  rsn ] = rotD100_func( eqe_selected_scaled_x[ rsn ].tolist(), eqe_selected_scaled_y[ rsn ].tolist())
        
        rotd100_mean_scaled_df["Mean"] = rotd100_mean_scaled_df[ rsn_selected ].mean( axis = 1 )
            
        plt.figure( figsize = (10 ,5))
        for rsn in key_list : 
            plt.plot( t , multiplied_selected_x[ rsn] , c = "gray" , lw = 0.5)
            plt.plot( t , multiplied_selected_y[ rsn] , c = "gray" , lw = 0.5)
        plt.plot( t , target["Sa"] , linestyle = "--" , c = "red", label= "Target")
        plt.plot( t , spectra_increament_factor *  target["Sa"] , c = "red" , lw = 1.5 ,  label = f"{spectra_increament_factor}x Target")
        plt.plot( t , geo_mean_1st_scaled_df[ "Mean"] , c = "blue" , label = "Geomean-Scaled")
        plt.plot( t , rotd100_mean_df[ "Mean"] , linestyle = "--" , c = "black" , label = "RotD100-Scaled")
        plt.plot( t , rotd100_mean_scaled_df[ "Mean"] , c = "black" , lw = 1.5 , label = "RotD100-Mean-Scaled")
        plt.legend( loc = 1)
        plt.axvspan( 0.2 * float( period) , 1.5 * float( period) , facecolor = "green" , alpha = 0.1)
        plt.xlim( 0  , 3 )
        plt.ylim(bottom = 0 )
        plt.box(False)
        plt.axhline(c="k")
        plt.axvline(c="k")
        plt.title("Geomean and RotD100 Spectra of the 2nd scaled selected ground motions")  ,
        plt.show()
        
        print( "-"*50)
        print( "Selected ground motions and scale factors")
        for rsn in rsn_selected : 
            temp_df = eqe_df[ eqe_df['RecordSequenceNumber'] == rsn]
            print( f"RSN{rsn} | { list( temp_df.EarthquakeName )[0] } | { sf_dict[ rsn ]}") 
        print( "-"*50)

    return(sf_dict)


def P695Scaling(key_list, target, period):

    # Selected Records DataFrame
    eqe_df = pd.read_csv("data/meta_data-R2.csv")

    spectral_data_x = pd.read_csv("data/spectral_x.csv")
    spectral_data_y = pd.read_csv("data/spectral_y.csv")

    selected_x = spectral_data_x.loc[spectral_data_x['RSN'].isin(key_list)]
    selected_y = spectral_data_y.loc[spectral_data_y['RSN'].isin(key_list)]

    t_str = selected_x.columns.tolist()[2:]
    t = []
    for i in t_str:
        t.append(float(i))

    rsn_selected = key_list

    eqe_selected_x = pd.DataFrame()
    eqe_selected_x.insert(0, 'T', t)
    eqe_selected_y = pd.DataFrame()
    eqe_selected_y.insert(0, 'T', t)

    for i in rsn_selected:
        eqe_selected_x[ i ] = selected_x.loc[selected_x['RSN'] == i].iloc[0].tolist()[2:]
        eqe_selected_y[ i ] = selected_y.loc[selected_y['RSN'] == i].iloc[0].tolist()[2:]

    ## Normalization Process ##

    # Create a DataFrame with PGV Values
    pgvFrame = pd.DataFrame(columns=['RSN', 'PGV'])
    pgvFrame['RSN'] = key_list

    pgvList = []
    try:
        for i in key_list:
            pgvList.append(eqe_df.loc[eqe_df['RecordSequenceNumber'] == i, 'PGV'].values[0])
    except:
        eqe_df = eqe_df.replace(';','', regex=True)
        eqe_df = eqe_df.rename(columns={'PGV;;;;': 'PGV'})
        for i in key_list:
            pgvList.append(eqe_df.loc[eqe_df['RecordSequenceNumber'] == i, 'PGV'].values[0])

    pgvFrame['PGV'] = pgvList

    pgvMedian = float(median(pgvList))

    normalizationFactors = []

    for i in pgvList:
        normalizationFactors.append(pgvMedian/float(i))

    pgvFrame['NM'] = normalizationFactors

    normalized_eqe_x = pd.DataFrame(columns=eqe_selected_x.columns)
    normalized_eqe_y = pd.DataFrame(columns=eqe_selected_y.columns)

    normalized_eqe_x['T'] = t
    normalized_eqe_y['T'] = t


    for i,j in zip(normalizationFactors, key_list):
        normalized_eqe_x[j] = [x * i for x in eqe_selected_x[j].to_list()]
        normalized_eqe_y[j] = [x * i for x in eqe_selected_y[j].to_list()]

    # Geometric Mean of Two Orientation # 
    geoMeanFrame = pd.DataFrame(columns=eqe_selected_x.columns)
    geoMeanFrame['T'] = t

    def geomean_func(acc_1 , acc_2):
        geo_mean = []
        for i , j in zip(acc_1 , acc_2):
            geo_mean.append(round((i * j)**(0.5), 4))
        return(geo_mean)

    for i in key_list:
        geoMeanFrame[i] = geomean_func(normalized_eqe_x[i], normalized_eqe_y[i])

    medianAcc = []

    for i in range(0, 111):
        medianAcc.append(median(geoMeanFrame.loc[i].values[1:].tolist()))

    geoMeanFrame['medianAcc'] = medianAcc


    # Median Spectral Acc. Value at Certain Period
    fMedian = interp1d(t, medianAcc)
    medianSa = fMedian(period)

    # Corresponding Target Spectral Acc.
    fTarget = interp1d(target['T'].to_list(), target['Sa'].to_list())
    targetSa = fTarget(period)

    # Determine Scale Factor
    scaleFactor = round((targetSa/medianSa), 3)

    plt.figure(figsize=(10,5))
    color_list = ("blue", "orange", "green", "red", "purple", "brown", "pink", "gray", "olive", "cyan", "lime")
    for i, c in zip(key_list, color_list):
        plt.plot(t, geoMeanFrame[i], linestyle = "--" , lw = .7 , color = c, label = i)
    plt.plot(t, target['Sa'], color='red', label='Target Spectra')
    plt.plot(t, geoMeanFrame['medianAcc'], color='blue', label = 'Median Spectra')
    plt.xlim(0  , 10)
    plt.ylim(bottom = 0 )
    plt.box(False)
    plt.axhline(c="k")
    plt.axvline(c="k")
    plt.axvline(x=period, linestyle = "--", label = 'Ta = {} sec'.format(period))
    plt.legend(loc = 1)
    plt.title("Median Spectra of Selected Set")
    plt.show()

    # Scaled Dataframe
    scaledFrame = pd.DataFrame(columns=eqe_selected_x.columns)
    scaledFrame['T'] = t

    for i in key_list:
        scaledFrame[i] = [x * scaleFactor for x in geoMeanFrame[i]]
    scaledFrame['medianAcc'] = [y * scaleFactor for y in geoMeanFrame['medianAcc']]

    plt.figure(figsize=(10,5))
    for i in key_list:
        plt.plot(t, scaledFrame[i], linestyle = "--" , lw = .7 , color = "gray")
    plt.plot(t, target['Sa'], color='red', label='Target Spectra')
    plt.plot(t, scaledFrame['medianAcc'], color='blue', label = 'Median Spectra Scaled')
    plt.xlim(0  , 10)
    plt.ylim(bottom = 0 )
    plt.box(False)
    plt.axhline(c="k")
    plt.axvline(c="k")
    plt.axvline(x=period, linestyle = "--", label = 'Ta = {} sec'.format(period))
    plt.legend()
    plt.title("Scaled Records")
    plt.show()

    print('Scale factor for selected data set is {}'.format(scaleFactor))