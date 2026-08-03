function outputSize = padSize(inputSize, minimumLength)
%PADSIZE Add trailing singleton dimensions if necessary.

arguments
    inputSize
    minimumLength (1,1) double {mustBeInteger, mustBePositive}
end

outputSize = inputSize;

if numel(outputSize) < minimumLength
    outputSize(end + 1:minimumLength) = 1;
end

end
