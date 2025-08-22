function D = readraw(fn, scanner)

if strcmp(scanner, 'GE')
    if strcmp(fn(end-2:end), '.h5')
        % ScanArchive file
        D = toppe.utils.loadsafile(fn, 'acq_order', true);
    else
        % P-file. Note flip along 1st (FID) dimension
        D = toppe.utils.loadpfile(fn, [], [], [], 'acqOrder', true, 'returnAsDouble', false); 
        D = flipdim(D,1);
    end
else
    twix = mapVBVD(fn);
    D = twix{2}.image.unsorted();
end
