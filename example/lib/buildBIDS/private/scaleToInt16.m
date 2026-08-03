function data = scaleToInt16(data, targetMagnitude)
%SCALETOINT16 Globally scale an array and convert it to int16.

arguments
    data
    targetMagnitude (1,1) double {mustBePositive}
end

data = double(data);

maximumMagnitude = max(abs(data(:)));

if maximumMagnitude > 0
    data = data .* (targetMagnitude / maximumMagnitude);
else
    warning('scaleToInt16:AllZeros', ...
        'Image contains only zeros; no scaling was applied.');
end

data = round(data);
data = min(data, double(intmax('int16')));
data = max(data, double(intmin('int16')));
data = int16(data);

end
