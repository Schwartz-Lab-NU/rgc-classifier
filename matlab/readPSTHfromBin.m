function data = readPSTHfromBin(fname, nfolds, loadData)

if nargin < 3
    loadData = true;
end

data = struct();

for f=1:nfolds
    fid = fopen(sprintf(fname,f),'rb');
    
    %save number of cells as int
    data(f).N = fread(fid,1,'uint16');
    data(f).nt = fread(fid,1,'uint16');
    data(f).ns = fread(fid,1,'uint16');
    
    
    %start by saving the start time, end time, and dt as double
    data(f).t0 = fread(fid,1,'double');
    data(f).tf = fread(fid,1,'double');
    
    data(f).dt = fread(fid,1,'double');
    
    %save smallest and largest spot sizes as int
    data(f).s0 = fread(fid,1,'uint16');
    data(f).sf = fread(fid,1,'uint16');
    
    %for each cell
%     d2 = d(folds{f});
%     l2 = l(folds{f});
%     minT2 = minT(folds{f});
%     maxT2 = maxT(folds{f});
%     se2 = se(folds{f});
%     allD2 = allD(folds{f},:,:);
    cs = struct();

    for n=1:data(f).N
        
        %save length of cellName as int, then cellName as char
        nchar = fread(fid,1,'uint8');
        cs(n).cellName = char(fread(fid, nchar, 'char'))';
        
        %save cellType as char (1 to 32)
        cs(n).label = fread(fid,1,'char');
        
        %save quality for each time point as logical
        cs(n).qt = fread(fid, data(f).nt,'logical');
        cs(n).qs = fread(fid, data(f).ns,'uint32');
        cs(n).psth = fread(fid,[data(f).nt, data(f).ns],'double');
        
        if ~loadData
            cs(n).psth = []; %save memory...
        end
    end
    
    data(f).cells = cs;
    
    fclose(fid);
end