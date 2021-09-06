import scalepy

record_df = scalepy.select_records('5 6', '0 250', '0 50', 'Strike - Slip')

#file_list = record_df['File Name (Horizontal 1)'].to_list()

#target_spectra = scalepy.target_spectrum(0.6, 0.4, "ZC")

#scalepy.SingleComponentScale(file_list, target_spectra, 2)