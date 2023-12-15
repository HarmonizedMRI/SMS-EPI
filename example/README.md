# SMS-EPI example from start to finish: acquisition and reconstruction

[under construction]

This folder contains a complete workflow for acquiring and reconstructing
SMS-EPI data for functional MRI in a vendor-agnostic manner.

The workflow involves the following steps:

1. **Create the Pulseq (.seq) files**
This example creates the following .seq files:
  * `cal.seq`: EPI ghost calibration scan
  * `acs.seq`: individual slice acquisitions for slice GRAPPA calibration
  * `sms-epi.seq`: SMS-EPI scan for fMRI
To create these sequence files, do:
```
>> main1;
```

TODO: blip up/down spin-echo scan for distortion correction

2. **Run the Pulseq files on your scanner**
At present, Siemens and GE scanners are supported.

3. **Reconstruct time-series images**
```
>> main2;
```


