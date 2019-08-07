# cell-migration-analysis
The driver code for this software is found in master.m

The driver code parses a spreadsheet with experimental data and does the 
following:

1) extractmetadata: extracts the experimental metadata from the ND2 file

2) extractRawData_alt: extracts the raw data from the tracking data

3) extractData: extracts data from rawdata and metadata

4) extractFeatures: extracts features (speed, persistence and angle) from the data

5) visualizeFeatures: displays various features chosen by the user

Note: to run the latest version (beta), run "master". To run the stable, less efficiently written version, run "master_version1"