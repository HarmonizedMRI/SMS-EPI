# SMS/3D EPI fMRI sequence for GE

Create 3D/SMS EPI fMRI sequence using the TOPPE library.


## Create the sequence files

```
>> [~,sys] = caipiepifmri('SMS', N, FOV, flip, nFrames, 'mb', 6);  % or '3D'
```

This creates the file 'fmri.tar' (in the Matlab working directory), that contains the following scan files:
```
  toppeN.entry
  seqstamp.txt
  scanloop.txt
  modules.txt
  .mod files
```
These files are described here: https://github.com/toppeMRI/toppe/blob/main/Files.md

These files are also written to the current Matlab working directory,
which allows you to plot the scan simply by typing:
```
>> toppe.plotseq(1,8,sys.ge);                % plot first 8 module executions
>> toppe.playseq(4,sys.ge, 'tpause', 0.2);   % loop through the entire scan
```

![Example sequence timing diagram](seq.png)

These files can also be **converted to Pulseq format** 
using the function `ge2seq.m` in the [https://github.com/toppeMRI/PulseGEq](PulseGEq) toolbox (WIP),
for execution on Siemens scanners.



## Reconstruct images

![Under (re)construction!](underreconstruction.jpg)

Scripts for loading data and gridding to cartesian grid are in utils/EPI.

<!--
Miscellaneous info

To create small jpg file (Linux):
cp ~/github/HarmonizedMRI/resourse/
inkscape -C -o underreconstruction.png ~/github/HarmonizedMRI/resource/images/underreconstruction.svg 
convert -quality 30 underreconstruction.png underreconstruction.jpg

-->
