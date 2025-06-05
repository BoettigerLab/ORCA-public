function AppendToDax(daxName,movie)

verbose = false;
binaryFormat = 'l';

fid = fopen(daxName, 'a');
if fid<0
    error(['Unable to open ' daxName]);
end

fwrite(fid, ipermute(uint16(movie), [2 1 3]), 'uint16', binaryFormat);

if verbose
    disp(['Finished appending ' daxName]);
end

fclose(fid);