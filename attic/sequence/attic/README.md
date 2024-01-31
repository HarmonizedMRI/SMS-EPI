# Cross-vendor SMS EPI fMRI sequence


## Dependencies and setup 

See [Setup.md](Setup.md)


## Create sequences

```
>> gre3d;      % for sensitivity maps (and B0 field mapping)
>> fmri2depi;  % (optional) for calibration and comparison
>> makesms;
```

## Create Pulseq and TOPPE scans 

```
>> makesms;
```
Creates `smsepi.seq` and `smsepi.tar`.

GE users:
Untar `SMSEPIfMRI.tgz` in /usr/g/bin/ and scan with the `toppev4` interpreter.

Siemens users:
Place `SMSEPIfMRI.seq` on scanner and execute scan.


