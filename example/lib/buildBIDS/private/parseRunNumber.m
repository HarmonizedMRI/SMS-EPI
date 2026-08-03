function runNumber = parseRunNumber(filename)
%PARSERUNNUMBER Extract the run number from task_runN.h5.nii.

arguments
    filename (1,:) char
end

tokens = regexp( ...
    filename, ...
    '^task_run(\d+)\.h5\.nii$', ...
    'tokens', ...
    'once');

if isempty(tokens)
    runNumber = [];
    return
end

runNumber = str2double(tokens{1});

if isnan(runNumber)
    runNumber = [];
end

end
