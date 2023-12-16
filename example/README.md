# A complete SMS-EPI example: from acquisition to image reconstruction

This folder contains a complete vendor-agnostic workflow for acquiring and reconstructing
SMS-EPI data for functional MRI.

The workflow involves the following steps:

1. **Create the Pulseq (.seq) files**.  
   ```
   >> main1;
   ```
   This creates the following .seq files:
   * `cal.seq`: EPI ghost calibration scan
   * `acs.seq`: individual slice acquisitions for slice GRAPPA calibration
   * `sms-epi.seq`: SMS-EPI scan for fMRI

   TODO: blip up/down spin-echo scan for distortion correction

2. **Run the Pulseq files on your scanner**.  
At present, Siemens and GE scanners are supported.

3. **Reconstruct time-series images**
   ```
   >> main2;
   ```

## TODO  
   * Siemens example data
   * Share data on DeepBlue/other
   * Can use similarity b/w calibration scans and recon'd images as QC metric


