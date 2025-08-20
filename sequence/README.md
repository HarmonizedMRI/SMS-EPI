This is github/jfnielsen/scanLog/EPI-test/nexttry/README.md

Create .seq files for SMS-EPI fMRI

.seq and .tar (TOPPE) scan files are located in 
UMich fMRI laboratory Google drive,
ENGIN-fMRI-laboratory > HarmonizedMRI project > fMRI > Protocol > sequene files > ABCD SMS-EPI

Oct 27, 2024: combine epi calibration and slice-grappa (mb=1) calibration scans into one
 ==> by simply setting 'doRefScan = true' for mb=1 scan. 
 This gives slice-by-slice ghost calibration data.

From early 2024:
Other recent changes:
  * getsmspulse.m: put first slice at -mb/2 (remove shift of +sliceSep/2)
     * writeEPI.m: accordingly, remove shift of -sliceSep/2 (rf.freqOffset)
  * set maxView = np x etl when mb>1, etl otherwise
  * Fixed slice offset for mb=1
  * etl=72 (multiple of mb=6)
  * TR=800ms
  * Add RF spoiling
  * add fat sat as default
  * Interleaved partition ordering, with last two even shots swapped
  * Remove arg.spoilersOn option (always on)

## Set up Python code for generating CAIPI sampling pattern

The ky-kz sampling pattern is generated using Rüdiger Stirnberg's Python code.
A copy of it is in

1. Get the code
Download from 
https://github.com/HarmonizedMRI/3DEPI.
On the Linux command line, do:
```
$ git clone git@github.com:HarmonizedMRI/3DEPI.git
```

2. Set up python virtual environment (recommended) and install required libraries in this environment:
```
sudo apt install python3.13-venv
python3 -m venv myvenv
source myvenv/bin/activate
pip install matplotlib
pip install scipy
```
