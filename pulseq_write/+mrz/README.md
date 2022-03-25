# function makeMBPulse
assumed conda is taking care of the python environment, otherwise remove this command, for example start MATLAB in the desired python environment.
calling example, band_sep is the center-to-center distance in terms of number of slices
> rf.signal = mrz.makeMBPulse(rf.signal,'n_bands',3,'timeBwProduct',4,'band_sep',11,'phs_0_pt','phs_mod');

# function makeMBPulse_less
no need to call python toolbox, **faster**
has more features than makeMBPulse
calling example
> rf.signal = mrz.makeMBPulse_less(rf.signal,'n_bands',0,'timeBwProduct',4,'band_sep',[0 -11 11]*5e-3/sliceThickness,'phs',[0 0.73 4.602]);%phs are correspondingly matched with slices
> rf.signal = mrz.makeMBPulse_less(rf.signal,'n_bands',2,'timeBwProduct',4,'band_sep',27*3e-3/sliceThickness);

