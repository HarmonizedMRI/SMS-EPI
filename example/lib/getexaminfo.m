function E = getexaminfo(fn, n)
% function E = getexaminfo(fn, n)
%
%  Inputs
%    fn     file name, i.e., 'sessions.txt'
%    n      row number

fid = fopen(fn, 'r');

% get file size
fseek(fid, 0, 1);
fsize = ftell(fid);
fseek(fid, 0, -1);

% skip to requested row
for i = 1:n-1
    l = fgetl(fid);
    assert(ftell(fid) < fsize, "end of file reached");
end

% read row
E = struct();
E.subject = fscanf(fid, '%s', 1);
E.site = fscanf(fid, '%s', 1);
E.scanner = fscanf(fid, '%s', 1);
E.date = fscanf(fid, '%s', 1);
E.session = fscanf(fid, '%s', 1);
E.vendor = fscanf(fid, '%s', 1);
E.kspace_delay = fscanf(fid, '%f', 1);
E.readout_trajectory_file = fscanf(fid, '%s', 1);

fclose(fid);

