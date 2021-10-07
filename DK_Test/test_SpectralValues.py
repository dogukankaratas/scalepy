import pandas as pd , os , numpy as np
from aad_TBDY2018_ResponseSpectrumCreator_class import tbdy2018_spectra

def test_spectral_values_1():

    spectrum_object = tbdy2018_spectra()
    spectrum_object.getSpectraValue(40.78823, 29.415702, "DD1")

    assert round(spectrum_object.spectral_value_dict['Ss'], 2) == 2.43
    assert round(spectrum_object.spectral_value_dict['S1'], 2) == 0.69
    assert round(spectrum_object.spectral_value_dict['PGA'], 2) == 0.96

def test_spectral_values_2():

    spectrum_object = tbdy2018_spectra()
    spectrum_object.getSpectraValue(38.240038, 26.806247, "DD1")

    assert round(spectrum_object.spectral_value_dict['Ss'], 2) == 2.12
    assert round(spectrum_object.spectral_value_dict['S1'], 2) == 0.54
    assert round(spectrum_object.spectral_value_dict['PGA'], 2) == 0.86

def test_spectral_values_3():

    spectrum_object = tbdy2018_spectra()
    spectrum_object.getSpectraValue(40.191556, 33.109492, "DD2")

    assert round(spectrum_object.spectral_value_dict['Ss'], 2) == 0.51
    assert round(spectrum_object.spectral_value_dict['S1'], 2) == 0.16
    assert round(spectrum_object.spectral_value_dict['PGA'], 2) == 0.22

def test_spectral_values_4():

    spectrum_object = tbdy2018_spectra()
    spectrum_object.getSpectraValue(39.74653, 39.521684, "DD2")

    assert round(spectrum_object.spectral_value_dict['Ss'], 2) == 1.47
    assert round(spectrum_object.spectral_value_dict['S1'], 2) == 0.42
    assert round(spectrum_object.spectral_value_dict['PGA'], 2) == 0.61

def test_spectral_values_5():

    spectrum_object = tbdy2018_spectra()
    spectrum_object.getSpectraValue(41.074643, 28.246586, "DD3")

    assert round(spectrum_object.spectral_value_dict['Ss'], 2) == 0.31
    assert round(spectrum_object.spectral_value_dict['S1'], 2) == 0.09
    assert round(spectrum_object.spectral_value_dict['PGA'], 2) == 0.14

def test_spectral_values_6():

    spectrum_object = tbdy2018_spectra()
    spectrum_object.getSpectraValue(40.5816068, 42.493853, "DD4")

    assert round(spectrum_object.spectral_value_dict['Ss'], 2) == 0.18
    assert round(spectrum_object.spectral_value_dict['S1'], 2) == 0.05
    assert round(spectrum_object.spectral_value_dict['PGA'], 2) == 0.08

