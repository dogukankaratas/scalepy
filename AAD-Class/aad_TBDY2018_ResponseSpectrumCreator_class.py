#%% IMPORT MODULES
# ımporting modules
# 
from ast import Break
from operator import index
import os, pandas as pd, matplotlib.pyplot as plt , numpy as np , datetime as dt , math  ,  scipy.interpolate as interp

from scipy.interpolate import interp1d


class tbdy2018_spectra:
    """
    This class estimates the spectral values of a given coordinate in accordance to Turkish Building Earthquake Code 2018. 
    The class derives the spectral values for TR grid from "AFAD_TDTH_parametre.csv" file. Please make sure you the file is in the same folder with this class file. 
    """

    # AFAD SPECTRAL VALUES
    """
    AFAD Turkish Seismic Hazard Spectral Values
    """
    try: 
        afad_spectra_params_df = pd.read_csv("AFAD_TDTH_parametre.csv")
    except Exception as err :
        print( f"{'*' * 20 }\n! Attention, AFAD  File is missing !\n{err}\n{'*' * 20 }")

    # PEER 
    """
    PEER FILES
    """
    PEER_FlatFile_path_name = "PEER_flat_final.xlsx"
    Spectra_file_x  = "PEER_records_spectra_x.csv"
    Spectra_file_y  = "PEER_records_spectra_y.csv"

    # __init__ 
    def __init__(self) :
        """
        These are the object attributes
        """
        self.spectral_value_dict = {}
        self.lat = 0
        self.lon = 0
        self.soil_class = ""
        self.intensity = ""
        self.period_list  = []
        self.spectral_orbits = []
        self.target_spectrum = []
        self.x_records_filtered_short = []
        self.y_records_filtered_short = []

    
    def getSpectraValue(self , lat , lon , intensity = "DD2"):
        """
        Interpolates the "Ss","S1","PGA","PGV" values from the "AFAD_TDTH_parametre.csv" file
        -Input-------------------------------------
        lat : latitude
        lon : longitude
        intensity : Earthquake Level; any of the given list "DD1" , "DD2" , "DD3" , "DD4". Defatuls is "DD2"

        -Output-------------------------------------
        spectral_value_dict : spectral values of "Ss","S1","PGA","PGV" at the given location. 
        """
        # intensity 
        self.intensity = intensity
        self.lat = lat 
        self.lon = lon

        # Grid locattions
        x = self.afad_spectra_params_df["LAT"].to_list()
        y = self.afad_spectra_params_df["LON"].to_list()
        
        # Spectral values dictionary
        for column_name in ["Ss","S1","PGA","PGV"]:

            z = self.afad_spectra_params_df[ f"{column_name}-{intensity}"].to_list()

            interpolator = interp.CloughTocher2DInterpolator( np.array([x,y]).T , z)

            spectral_value = np.round( interpolator( lat,lon)  , 3 )
            
            self.spectral_value_dict[ column_name ] = spectral_value

    def show_spectral_values(self):
        """
        Printing the spectral values estimated at the given location.         
        """
        [ print( f"{column_name}-{self.intensity} { item }") for column_name , item in self.spectral_value_dict.items() ]


    def get_spectral_ordinates( self , vs30 ) :
        """
        -Input-------------------------------------
        vs30 : Soil Velocity in m/s units
        -Output-------------------------------------
        self.period_list : Period values 
        self.orbital_values : RS values
        """

        # Determine the short and 1 sec spectral values. 
        if self.spectral_value_dict == {} : 
            print( "Provide location coordinates to estimate the spectral values, first. Then run self.get_spectral_ordinates() method.")
            exit()
        else:
            Ss = self.spectral_value_dict["Ss"] 
            S1 =  self.spectral_value_dict["S1"] 
        # zemin bilgisi belirlensin
        #soil_class = self.getSoilClass( vs30 )

        vs30_values = [ 0 , 180 , 360 , 760 , 1_500 , 20_000 ]
        soil_class_list = [ "ZE" , "ZD" , "ZC" , "ZB" , "ZA" ]
        vs_limit , count = 0 , 0
        while vs30 >= vs_limit:
            soilClass =  soil_class_list[ count ]
            count += 1
            vs_limit = vs30_values[ count]
        # Class property
        self.soil_class = soilClass


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

        # Short period estimation
        if Ss < Ss_range[0]:
            Fs = FS_table[self.soil_class][0]
            SDs = Ss * Fs
        elif Ss > Ss_range[-1]:
            Fs = FS_table[self.soil_class][-1]
            SDs = Ss * Fs    
        else:
            FS_satir = interp1d(Ss_range, FS_table[self.soil_class], kind='linear')
            FS_katsayisi = FS_satir(Ss)
            Fs = round( float(FS_katsayisi) , 2) 
            SDs = Ss * Fs
        # 1sn period estimation
        if S1 < S1_range[0] :
            F1 = F1_table[self.soil_class][0]
            SD1 = S1 * F1
        elif S1 > S1_range[-1]:
            F1 = F1 = F1_table[self.soil_class][-1]
            SD1 = S1 * F1
        else:    
            F1_satir = interp1d(S1_range, F1_table[self.soil_class], kind='linear')
            F1_katsayisi = F1_satir(S1)
            F1 = round(float(F1_katsayisi) , 2)
            SD1 = S1 * F1

        # Corner period values
        TA = 0.2 * SD1 / SDs
        TB = SD1 / SDs
        TL = 6

        # Horizontal acceleration ordinates
        def spektra_yatay(T,SDs,SD1, TA, TB , TL):  
            if T < TA :
                return((0.4 + 0.6*(T/TA))*SDs)
            elif T >= TA and T <= TB:
                return(SDs)
            elif T> TB and T <= TL:
                return(SD1 / T)
            elif T> TL:
                return(SD1*TL/(T**2))

        # Spectrum values
        self.period_list  = [ item  for item in np.arange(0, 6.05, 0.076)]

        self.spectral_orbits = [ round( spektra_yatay(period,SDs,SD1, TA, TB , TL) , 3 ) for period in self.period_list ]


    def spectra_plot( self , plot_save = False):
        """
        Visualization of spectra
        """
        if self.spectral_orbits == [] :
            print( "Please run self.get_spectral_ordinates() method first")
            raise()
        else: 
            plt.figure()
            plt.axhline(c="black")
            plt.axvline(c="black")        
            plt.plot( self.period_list , self.spectral_orbits)
            plt.xlabel( "Period (s)")
            plt.ylabel("Sa (g)")
            plt.title( f"{ self.intensity} - {self.soil_class}")
            plt.box( False)
            plt.xlim( 0 , self.period_list[-1])
            plt.ylim(bottom = 0 )
            plt.tight_layout()
            if plot_save == True : 
                plt.savefig( f"{self.intensity}-{self.soil_class}-{self.lat}_{self.lon}-RS.png")
                plt.show()
            else :
                plt.show()

    def select_records_aad( self, magnitude_range, vs_range, rjb, fault_mechanism, T1):

        period_list = [ item  for item in np.arange(0, 6.05, 0.076)]

        eqe_df = pd.read_excel(self.PEER_FlatFile_path_name )

        ## Reading the ranges
        min_m, max_m = [float(x) for x in magnitude_range.replace(" ","").split(",")]
        min_vs, max_vs = [float(x) for x in vs_range.replace(" ","").split(",")]
        min_d, max_d = [float(x) for x in rjb.replace(" ","").split(",")]

        ## Filtering the earthquakes
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


        # Removing not used data ##
        eqe_selected.drop(eqe_selected.columns.difference(['Record Sequence Number', 'Earthquake Name', 
                                                        'File Name (Horizontal 1)', 'File Name (Horizontal 2)']), 1, inplace = True)
        # Column renaming #                                             
        eqe_selected = eqe_selected.rename(columns={'Record Sequence Number': 'ID'})

        # Cleaning #
        del eqe_s
        del eqe_df

        # Reading PEER files #
        x_records = pd.read_csv( self.Spectra_file_x , index_col = False)
        y_records = pd.read_csv( self.Spectra_file_y , index_col = False)


        # Spectral values assigned #
        x_records_filtered = x_records.loc[:, x_records.columns.isin(eqe_selected['ID'].astype(str).to_list())]
        y_records_filtered = y_records.loc[:, y_records.columns.isin(eqe_selected['ID'].astype(str).to_list())]

        # Geomean of the components #
        geo_mean_df = pd.DataFrame()
        geo_mean_df["T"]  = period_list

        for i in x_records_filtered.columns:
            geo_mean_df[i] = [(x*y)**(1/2) for x,y in zip(x_records_filtered[str(i)].to_list(), y_records_filtered[str(i)].to_list())]

        # For the given period range estimate the Scale Factor by minimizing the difference
        geo_mean_range_df = geo_mean_df[ ( geo_mean_df["T"] >= 0.2 * T1) & ( geo_mean_df["T"] <= 1.5 * T1) ]
        
        data = {'T': self.period_list  , 'Sa': self.spectral_orbits }  
        self.target_spectrum = pd.DataFrame().from_dict(data ) 

        target_spectrum_range_df = self.target_spectrum[ ( self.target_spectrum["T"] >= 0.2 * T1) & ( self.target_spectrum["T"] <= 1.5 * T1) ]

        ## Selection of the 11 records ##
        difference_dict  = {}

        for column_name in geo_mean_df.columns[1:]:
            difference_dict[ column_name ] = sum( geo_mean_range_df[ column_name ] * target_spectrum_range_df["Sa"] ) / sum(target_spectrum_range_df["Sa"]**2 )


        eqe_selected["Diff_in_range"] = difference_dict.values()

        eqe_selected_sorted = eqe_selected.sort_values(by ="Diff_in_range" )

        first_three_list = []

        for counter, temp_df in eqe_selected_sorted.groupby( by = "Earthquake Name") :
            [first_three_list.append(item) for item in  temp_df["ID"][:3].to_list() ]

        eqe_selected_sorted_short = eqe_selected_sorted[ eqe_selected_sorted["ID"].isin( first_three_list )]
        eqe_selected_sorted_short.sort_values(by = "Diff_in_range", inplace= True)

        eqe_selected_sorted_short = eqe_selected_sorted_short.iloc[:11 , : ]

        if len( eqe_selected_sorted_short ) != 11: 
            print( f"{'x'*50}\n! Not enough earthquake record available !\n{'x'*50}")
            exit( )

        # Spectral Values 

        self.x_records_filtered_short = x_records_filtered[ [str( item) for item in eqe_selected_sorted_short["ID"].to_list()]]
        self.y_records_filtered_short = y_records_filtered[ [str( item) for item in eqe_selected_sorted_short["ID"].to_list()]]
        #return( eqe_selected_sorted_short , self.x_records_filtered_short , self.y_records_filtered_short)

        
    #def BothComponentScaling(self, records_df_x , records_df_y, target_spectra , T1):   
    def BothComponentScaling(self , T1):   

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

        spectra_df = pd.DataFrame()

        T = np.arange(0.05, 4.05, 0.05)
        spectra_df['Periyot'] = T

        for i, j in zip( self.x_records_filtered_short.columns.to_list(), self.y_records_filtered_short.columns.to_list()):
            spectra_df["RSN" + str(i) + "_000"] = self.x_records_filtered_short[str(i)]
            spectra_df["RSN" + str(j) + "_090"] = self.y_records_filtered_short[str(j)]
        
        spectra_geomean_df = pd.DataFrame()
        spectra_geomean_df["Periyot"] = spectra_df["Periyot"]
        kolon_isimleri = spectra_df.columns[1:]
        for i,j in zip( kolon_isimleri[::2] , kolon_isimleri[1::2]):
            #print( i , j )
            nameEQE = i.split("_")[0]
            spectra_geomean_df[nameEQE] = geomean_func( spectra_df[i] , spectra_df[j])

        carpim = {}
        for column in spectra_geomean_df.columns[1:] :
            plt.plot( spectra_geomean_df["Periyot"] , spectra_geomean_df[ column] ,  "gray", lw= 0.5 )
            carpim[ column ] = round( sum( ( self.target_spectrum["Sa"]  * spectra_geomean_df[ column]) ) / sum( spectra_geomean_df[ column]**2 ) , 3)

        spectra_srss_mean_df = pd.DataFrame()
        spectra_srss_mean_df["Periyot"] = spectra_df["Periyot"]
        kolon_isimleri = spectra_df.columns[1:]

        plt.figure( figsize= [ 10 , 5 ])

        for i,j,carpan in zip( kolon_isimleri[::2] , kolon_isimleri[1::2] , carpim.keys() ) :
           # print( i , j ,carpan , carpim[carpan] )
            nameEQE = i.split("_")[0]
            spectra_srss_mean_df[nameEQE] = srss_func( carpim[carpan]*spectra_df[i] , carpim[carpan]*spectra_df[j] )
            plt.plot( spectra_srss_mean_df["Periyot"] ,  spectra_srss_mean_df[nameEQE] , "gray", lw= 0.3)

        ortalama = spectra_srss_mean_df.iloc[:,1:].mean( axis = 1 )

       
        T1_alt , T1_ust = round( 0.2 * T1  ,2 ) , round(1.5*T1  , 2)

        [ plt.axvline( item  , color= "k", linestyle = "-." ) for item in [T1 , T1_alt , T1_ust] ]

        indexAralik = spectra_srss_mean_df.index[(spectra_df["Periyot"] >= T1_alt) & (spectra_srss_mean_df["Periyot"] <= T1_ust)].tolist()

        plt.axvspan( T1_alt , T1_ust  , alpha=0.1, color='blue')

        spectra_increament_factor = 1.3 

        # Scale Öncesi
        plt.plot( spectra_srss_mean_df["Periyot"] , spectra_geomean_df.iloc[:,1:],"c" , lw = 0.3 )

        plt.plot( spectra_srss_mean_df["Periyot"] , spectra_geomean_df.iloc[:,1:].mean(axis=1),"r--" , lw = 1 , label="Non-Scaled Average" )

        plt.plot( self.target_spectrum["T"] , self.target_spectrum["Sa"], "blue" , lw= "2" , label="Target Spectrum" )

        plt.plot( self.target_spectrum["T"] , spectra_increament_factor * self.target_spectrum["Sa"], "b--" , lw= "2" , label="1.3 x Target Spectrum" )

        plt.plot( self.target_spectrum["T"] , ortalama , "r-" ,lw= "1" , label="Preliminary Scale Average" )
        plt.legend()

        SF_ortalama , tol , counter  = 1 , 1 , 0

        while tol > 0.01:
            
            ortalama_Scaled = ortalama.iloc[indexAralik] * SF_ortalama 

            farklar = spectra_increament_factor * self.target_spectrum["Sa"].iloc[ indexAralik ] - ortalama_Scaled 
            
            if max(farklar) > 0 : 
                SF_ortalama  += 0.01
            
            if max(farklar) < 0 :
                SF_ortalama  -= 0.01
            
            if max( farklar ) > 0 and  max(farklar) < 0.01 : 
                tol = 0.001 

            counter += 1
            if counter == 50: 
                tol = 0.001

        plt.plot( self.target_spectrum["T"] , ortalama * SF_ortalama , "green" ,lw= "2" , label="Last Scaled Average" )

        plt.xlabel("Period (s)"),plt.ylabel("Sa (g)")
        plt.legend()
        #plt.xlim(left= 0 , right = self.target_spectrum["T"].iloc[-1]), plt.ylim( bottom = 0 )
        plt.xlim(left= 0 , right = 4 ), plt.ylim( bottom = 0 )

        plt.title( f"Scaled Spectra for T={T1}" )
        plt.tight_layout()

        #%% Ölçek Katsayısını hesaplanması
        ScaleFactors = {}
        for key in carpim.keys():
            ScaleFactors[ key ] = round( carpim[key] * SF_ortalama , 3) 

        ScaleFactors = dict( sorted( ScaleFactors.items() ) ) 

        for key, value in ScaleFactors.items():
            print( f"{key} SF = {value}")