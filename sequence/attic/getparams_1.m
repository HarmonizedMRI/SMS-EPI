function [seq, system] = getparams
%
% Define common acquisition parameters and system hardware parameters

seq.nx = 92; seq.ny = seq.nx;   % matrix size
seq.fov = 22.08;                 % in-plane fov (cm)
seq.slThick = seq.fov/seq.nx;   % isotropic voxels

seq.nSpoilCycles = 1.5;           % lower bound on number of cycles of gradient spoiling across voxel (applied to x and y)

% Scan plane prescription (rotation and slice offset)
if 0
	% use slab prescription from SlicePlanner GUI (https://github.com/toppeMRI/SlicePlanner)
	roi = toppe.getroi('ROI.h5', 1); 
	seq.rotmat = roi.rotmat;
	resLoc = 0.1; % Voxel size of localizer volume (assumed to be isotropic) (cm).            
	seq.sliceOffset = resLoc*roi.scanPlaneToIsocenterDistance;      % cm
else
	% scan at iso-center, non-oblique
	seq.rotmat = eye(3);
	seq.sliceOffset = 0;
end

% Hardware limits
% NB! maxGrad must be equal to physical hardware limit since used to scale gradients.
% maxSlew can be less than physical limit (for reducing PNS)
system.ge = toppe.systemspecs('maxSlew', 20, 'slewUnit', 'Gauss/cm/ms', ...
	'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...  
	'maxRF', 0.25, 'rfUnit', 'Gauss');
system.ge.forbiddenEspRange = [0.41 0.51];    % (ms) Forbidden echo spacings (mechanical resonance). See /usr/g/bin/epiesp.dat
system.ge.toppe.myrfdel = 78;        % rf/gradient delay (us) ( = psd_rf_wait). Inside scanner: 78. Outside scanner: 94.
system.ge.toppe.daqdel  = 84;        % daq/gradient delay (us) (= psd_grd_wait). Inside scanner: 84. Outside scanner: 100.

%system.siemens = mr.opts('MaxGrad', 28, 'GradUnit', 'mT/m', ...
%    'MaxSlew', 150, 'SlewUnit', 'T/m/s', ...
%	 'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

% Fat saturation pulse (bandwidth = 500 Hz)
seq.fatsat.flip    = 90;
seq.fatsat.slThick = 1000;       % dummy value (determines slice-select gradient, but we won't use it; just needs to be large to reduce dead time before+after rf pulse)
seq.fatsat.tbw     = 1.5;
seq.fatsat.dur     = 3;          % pulse duration (msec)
seq.fatsat.freqOffset = -440;    % Hz

