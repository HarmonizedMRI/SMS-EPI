# Cross-vendor SMS EPI fMRI sequence


## Dependencies and setup 

See [Setup.md](Setup.md)


## Test the code

```
>> cd SMS-EPI/sequence/
>> makeSMSEPIfMRI('test');
```

## Create Pulseq and TOPPE scans 

```
>> makeSMSEPIfMRI;
```
Creates `SMSEPIfMRI.seq` and `SMSEPIfMRI.tgz`.

GE users:
Untar `SMSEPIfMRI.tgz` in /usr/g/bin/ and scan with the `toppev4` interpreter.

Siemens users:
Place `SMSEPIfMRI.seq` on scanner and execute scan.


