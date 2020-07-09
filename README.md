# SiPM_waveform_analysis
****************************************************************************
Tools to process waveform in binary format acquired with LeCroy oscilloscope
****************************************************************************

0) Copy binary scope data to "data" folder



1) Convert scope data from binary to root using:

> python convertWaveToRoot.py folder_name_containing_WF  prefix_of_waveforms output_file_name

--> output raw root files are saved in the folder: "root_files"



2) Analyze waveforms to extract time stamps with different methods:

> ./pulses_analyzer.exe input_file_path output_file_name (bias optional)


--> output summary root files are saved in the folder: "output"


