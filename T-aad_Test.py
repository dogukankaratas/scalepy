"""
"""
#%%
import scalepy as sp

#%%
target = sp.targetSpectrum(0.8, 0.4, 'ZE')

record_keys = sp.recordSelection('6 8', '100 280', '0 250', 'Strike - Slip', '0 15', '0 20', '0 5', target , period = 1 ,)

RSNList_SF_dict = sp.amplitudeScaling(record_keys, target, period = 1 )

# %%
