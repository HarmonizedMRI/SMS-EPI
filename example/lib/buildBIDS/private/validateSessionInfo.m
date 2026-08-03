function validateSessionInfo(info)
%VALIDATESESSIONINFO Validate fields required by modality-copy functions.

arguments
    info (1,1) struct
end

requiredFields = {'srcDir', 'sub', 'ses', 'sessionDir'};

for iField = 1:numel(requiredFields)
    fieldName = requiredFields{iField};

    if ~isfield(info, fieldName) || isempty(info.(fieldName))
        error('validateSessionInfo:MissingField', ...
            'Session info is missing required field: %s', ...
            fieldName);
    end
end

end
