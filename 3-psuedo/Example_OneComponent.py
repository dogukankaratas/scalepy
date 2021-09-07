import scalepy

target_spectra = scalepy.target_spectrum(0.6, 0.4, "ZC")

record_df = scalepy.select_records('5 6', '0 250', '0 50', 'Strike - Slip')

