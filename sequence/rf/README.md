# Make SMS RF pulse, and sequences for imaging slice profile

## Dependencies and setup 

See [Setup.md](Setup.md)


## Test the code

```
>> cd SMS-EPI/sequence/rf/
>> makeSMSpulse('test');
```

## Create Pulseq and TOPPE scans for imaging the slice profile

```
>> makeProfileScan;
```
Creates `SMSprofile.seq` and `SMSprofile.tgz`.

GE users:
Untar `SMSprofile.tgz` in /usr/g/bin/ and scan with the `toppev3` interpreter.

Siemens users:
Place `SMSprofile.seq` on scanner and execute scan.


