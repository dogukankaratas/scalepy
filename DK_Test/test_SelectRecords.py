import pandas as pd , os , numpy as np
from aad_TBDY2018_ResponseSpectrumCreator_class import tbdy2018_spectra


vs30_list = [200, 240, 280, 320, 360]
magnitude_range = "4,9" 
vs_range = "180,360"
rjb_range = "50,100"
fault_list = [ "Strike - Slip" ,  "Normal" , "Reverse" , 'Reverse - Oblique' , 'Normal - Oblique']
T1 = .2

lat_list = [40.78823, 38.240038, 40.191556, 39.74653, 41.074643]
long_list = [29.415702, 26.806247, 33.109492, 39.521684, 28.246586]


for lat,lon in zip(lat_list, long_list):
        
    spectrum_object = tbdy2018_spectra()
    spectrum_object.getSpectraValue(lat, lon, "DD2")
    spectrum_object.get_spectral_ordinates(320)
    spectrum_object.select_records_aad( magnitude_range, vs_range, rjb_range, "Reverse", T1)
    spectrum_object.BothComponentScaling(T1)

        
