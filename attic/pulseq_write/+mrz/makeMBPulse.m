function [rf2] = makeMBPulse(pulsein,varargin)
%makeMBPulse make an SMS pulse
%     a wrapper to a python function(see below). See supported params below
%     in the 'parser' section. 
%     work on windows also, anaconda is convenient for customerizing the
%     environment, same with Linux. Install the required Python library by
%     executing "pip3 install sigpy" or "pip install sigpy", if not 'pip',
%     there are other ways installing the toolbox, refer to the document.
%     mac not tested yet, no toy around so far.
%
%     BE CAREFUL, MB stands for the method multi rf added together
%
%       mb_rf(pulse_in, n_bands=3, band_sep=5, phs_0_pt='None')
%          Args:
%              pulse_in (array): samples of single-band RF pulse.
%              n_bands (int): number of bands.
%              band_sep (float): normalized slice separation.
%              phs_0_pt (string): set of phases to use. Can be 'phs_mod' (Wong),
%                 'amp_mod' (Malik), 'quad_mod' (Grissom), or 'None'
% 
%          band_sep = slice_sep/slice_thick*tb, where tb is time-bandwidth product
%          of the single-band pulse
% 
%          Returns:
%              array: multibanded pulse out
% 
%          References:
%              Wong, E. (2012). 'Optimized Phase Schedules for Minimizing Peak RF
%              Power in Simultaneous Multi-Slice RF Excitation Pulses'. Proc. Intl.
%              Soc. Mag. Reson. Med., 20 p. 2209.
%              Malik, S. J., Price, A. N., and Hajnal, J. V. (2015). 'Optimized
%              Amplitude Modulated Multi-Band RF pulses'. Proc. Intl. Soc. Mag.
%              Reson. Med., 23 p. 2398.
%
%       PINS is not supported yet
%     Periodic modulation RF strategies, such as PINS pulses, have an un-
% limited FOV in the slice direction, requiring a careful selection of slice
% orientation, which may not necessarily be optimal for slice accelerations.
% Alternative approaches include variable-rate selective excitation
% (VERSE) (Conolly et al., 1988; Setsompop et al., 2012), optimized phase
% offsets between bands (Goelman, 1997; Hennig, 1992; Wong, 2012;
% Yao, 1995), or time shifting between individual bands can be applied
% (Auerbach et al., 2013; Goelman, 1997; Yao, 1995) to reduce peak power.
%       pins(tb, sl_sep, sl_thick, g_max, g_slew, dt, b1_max=0.18, ptype='ex', 
%       ftype='ls', d1=0.01, d2=0.01, gambar=4258)
%         Parameters:
%             • tb (float) – time-bandwidth product. 
%             • sl_sep (float) – slice separation in cm. 
%             • sl_thick (float) – slice thickness in cm. 
%             • g_max (float) – max gradient amplitude in gauss/cm 
%             • g_slew (float) – max gradient sliew in gauss/cm/s 
%             • dt (float) – RF + gradient dwell time in s. 
%             • b1_max (float) – Maximum RF amplitude 
%             • ptype (string) – pulse type, ‘st’ (small-tip excitation), ‘ex’ (pi/2 excitation pulse), ‘se’ (spin-echo pulse), ‘inv’ (inversion), or ‘sat’ (pi/2 saturation pulse). 
%             • ftype (string) – type of filter to use in pulse design 
%             • d1 (float) – passband ripple level in :math:’M_0^{-1}’. 
%             • d2 (float) – stopband ripple level in :math:’M_0^{-1}’. 
%             • gambar (float) – Appropriate gyromagnetic ratio in Hz/gauss. 
%         Returns:
%             • rf (array): RF Pulse in Gauss 
%             • g (array): Gradient waveform in Gauss/cm 
%         
%         References:
%             Norris, D.G. and Koopmans, P.J. and Boyacioglu, R and Barth, M (2011).
%             'Power independent of number of slices (PINS) radiofrequency Pulses
%             for low-power simultaneous multislice excitation'.
%             Magn. Reson. Med., 66(5):1234-1240.
%
%       2022.03.08

            
% validPulseTypes = {'pins','pha_offset'};
validPulseUses = mr.getSupportedRfUse();

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeMBPulse';
    
    % RF params
    addRequired(parser, 'pulsein', @isnumeric);
    addOptional(parser, 'system', mr.opts(), @isstruct);
    addParamValue(parser, 'timeBwProduct', 4, @isnumeric);   
    addParameter(parser,'phs_0_pt','None',@(x) any(validatestring(x,{'None','phs_mod','amp_mod','quad_mod'})));
    addParameter(parser,'n_bands',3,@isnumeric);
    addParameter(parser,'band_sep',5,@isnumeric);
    % Slice params
    addParamValue(parser, 'maxGrad', 0, @isnumeric);
    addParamValue(parser, 'maxSlew', 0, @isnumeric);
    addParamValue(parser, 'sliceThickness', 0, @isnumeric);
    addParamValue(parser, 'delay', 0, @isnumeric);
    addParamValue(parser, 'dwell', 0, @isnumeric); % dummy default value    
    % whether it is a refocusing pulse (for k-space calculation)
    addOptional(parser, 'use', '', @(x) any(validatestring(x,validPulseUses)));
end

parse(parser, pulsein, varargin{:});
opt = parser.Results;

if opt.dwell==0
    opt.dwell=opt.system.rfRasterTime;
end
nBands  = opt.n_bands;
tb      = opt.timeBwProduct;
BandSep = opt.band_sep;
method  = opt.phs_0_pt;


% mac not tested yet
% env is a bit trickier, starting matlab from an activated python env will
% do all the trick already, otherwise the env should be activated during
% the system calling.
% python = 'C:\Users\zhangs\Documents\anaconda3\python.exe';
% conda = 'C:\Users\zhangs\Documents\anaconda3\Scripts\activate';
if ispc
    [status, python] = system('where python');
else
    [status, python] = system('which python');
end
if ~status
    python = strsplit(python);
    python = python{1};
else
    error('please set the python path')
end
if ispc
    [status, conda] = system('where activate');
else
    [status, conda] = system('which activate');
end
if ~status
    conda = strsplit(conda);
    conda = conda{1};
else
    warning('conda is not found, env may not be activated and ready')
end

% %pins, discard since not efficient
% %       pins(tb, sl_sep, sl_thick, g_max, g_slew, dt, b1_max=0.18, ptype='ex', 
% %       ftype='ls', d1=0.01, d2=0.01, gambar=4258)
% code = ['import sigpy.mri.rf.multiband as mb\n' ...
%     'pulse= mb.dz_pins(tb=4,sl_sep=1.2,sl_thick=0.3,g_max=3,g_slew=15000,dt=1e-5,b1_max=0.18,ptype="st",ftype="ls",d1=0.01,d2=0.01,gambar=4258)\n' ...
%     'print(*pulse[0])\n' 'print(*pulse[1])'];
% [status,results] = system(cmd);
% rf3 = str2num(results);

tmp = num2str(pulsein);% env cannot be longer than 2^15, not help too much tmp1 = num2str(rf.signal,4);
strlen = ceil(length(tmp)/32766);
% for stri = 1:strlen
%     if stri == strlen
%         eval(['val' num2str(stri) ' = tmp((stri-1)*32766+1:end);'])
%     else
%         eval(['val' num2str(stri) ' = tmp((stri-1)*32766+1:stri*32766);'])
%     end
%     eval(['key' num2str(stri) ' = ''hao' num2str(stri) ''';'])
% end
val = cell(strlen,1);
key = cell(strlen,1);
for stri = 1:strlen
    if stri == strlen
        val{stri} = tmp((stri-1)*32766+1:end);
    else
        val{stri} = tmp((stri-1)*32766+1:stri*32766);
    end
    key{stri} = ['hao' num2str(stri)];
end
for stri = 1:strlen
    setenv(key{stri},val{stri})
end
for stri = 1:strlen
    if stri == 1
        env = ['os.getenv(''hao' num2str(stri) ''')'];
    else
        env = [env '+' 'os.getenv(''hao' num2str(stri) ''')'];
    end
end
% Linux way
% bash: python -c $'import os\nimport numpy as np\nnp.savetxt("tmp",[12,3])'
% bash: python -c "import os;import numpy as np; np.savetxt('tmp',[12,3])"
% Windows way
% batch: activate & python -c "exec(\"import numpy;\nimport matplotlib;\nnumpy.savetxt('tmp',[12,3])\")"
% batch: activate & python -c "import numpy;import matplotlib;numpy.savetxt('tmp1',[12,3]);"
cmd = [conda ' & ' python ' -c "import sigpy.mri.rf.multiband as mb;' ...
    'import numpy as np;import os;'...
    'pulsein = ' env ';'...
    'pulsein = pulsein.split();'...
    'pulsein = [complex(i) for i in pulsein];'...   if contains complex 
    ];
cmd2= [
    'rf = mb.mb_rf(pulsein,n_bands=%d,band_sep=%d,phs_0_pt=''%s''' ');' ...
    ...'np.savetxt(''tmp''' ',pulsein.split('' ''' '))'
    'print(*rf)"'...
    ];
[status,results] = system([cmd,sprintf(cmd2,nBands,tb*BandSep,method)]);
if status~=0
    display(results)
    error('executing python command failed');
end

rf2 = str2num(results);
for stri = 1:strlen  %to delete the env
    setenv(key{stri},'')
end


end%function end

