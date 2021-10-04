
#%%
import pandas as pd , os , numpy as np
from aad_TBDY2018_ResponseSpectrumCreator_class import tbdy2018_spectra


#%%
"""
fault_list = [ "Strike - Slip" ,  "Normal" , "Reverse" , 'Reverse - Oblique' , 'Normal - Oblique']
"""


vs30 = 250 
magnitude_range = "4,9" 
vs_range = "180,360"
rjb_range = "50,100"
fault_mechanism = "Reverse"
T1 = .2

spectrum_object = tbdy2018_spectra()
spectrum_object.getSpectraValue( 39.39062046395253, 32.307622830661735 , "DD2" )
spectrum_object.show_spectral_values()
spectrum_object.get_spectral_ordinates( vs30 )
spectrum_object.select_records_aad( magnitude_range, vs_range, rjb_range, fault_mechanism, T1)
spectrum_object.BothComponentScaling(T1)

# %%
