# SMS/3D EPI fMRI sequence

# Create scan files for GE

```
>> [seq,sys] = caipiepifmri('SMS');  % or '3D'
>> toppe.plotseq(1,8,sys.ge);
>> toppe.playseq(4,sys.ge, 'tpause', 0.2);
```
