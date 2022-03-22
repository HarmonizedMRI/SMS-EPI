function sys = getsys
% function sys = getsys
%
% Return 'sys' struct containing GE and Siemens hardware specs.
%
% Edit as needed for your scanner.

% Set hardware limits (for design and detailed timing calculations).
% 'maxSlew' and 'maxGrad' options can be < scanner limit, and can vary across .mod files. 

% system struct for GE
% maxSlew and maxGrad can be less than physical limits
sys.ge = toppe.systemspecs('maxSlew', 20, 'slewUnit', 'Gauss/cm/ms', ...
    'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
    'myrfdel', 200, ...  % psd_rf_wait
    'daqdel', 100, ...   % psd_grd_wait
    'timessi', 100, ...
    'gradient', 'xrm');  % MR750 scanner. Used to plot PNS profile in toppe.plotmod()

% system struct for Siemens
sys.siemens = mr.opts('MaxGrad', sys.ge.maxGrad*10, 'GradUnit', 'mT/m', ...
    'MaxSlew', sys.ge.maxSlew*10, 'SlewUnit', 'T/m/s', ...
    'B0', 3.0, ...
    'rfRingdownTime', 20e-6, ...
    'rfDeadTime', 100e-6, ...
    'adcDeadTime', 10e-6);

