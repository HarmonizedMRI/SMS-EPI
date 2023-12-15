# SMS-EPI example from start to finish: acquisition and reconstruction

[under construction]


## Overview

This folder contains a complete workflow for acquiring and reconstructing
SMS-EPI data for functional MRI.

The workflow involves the following steps:

1. **Create the Pulseq (.seq) file**
This is done by scanning a uniform phantom with a Pulseq sequence that we provide 
(see the [sequence/Pulseq](sequence/Pulseq) folder), 

2. **Run the Pulseq file on your scanner**
This involves running a Pulseq scan.
At present, Siemens and GE scanners are supported.

3. **Load the raw data**
We provide MATLAB scripts for reading Siemens .dat files,
and GE P-files or Raw Data Server (RDS) files.

4. **Reconstruct time-series images**


