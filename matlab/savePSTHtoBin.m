function [allT, allS] = savePSTHtoBin(d, dt, fname, folds, allT, allS)
%d is a struct array (e.g. table2struct(single_SMS_table))


%rebin the data to a common frame -- expansion
t = {d(:).PSTH_X};
minT = cellfun(@min,t);
maxT = cellfun(@max,t);

if (nargin < 5)
    allT = min(minT):dt:max(maxT);
end

s = {d(:).SpotSizeVec};

if (nargin < 5)
    allS = min(cellfun(@min,s)):max(cellfun(@max,s));
end

[gt,gs] = meshgrid(allT,allS);

allD = reshape(rebinData({d(:).SMS_PSTH},s,gs(:)',[minT' maxT'],gt(:),dt), [], numel(allS), numel(allT));

se = cellfun(@(x) min(abs(allS' - x),[],2).^2,s,'uniformoutput',false);
[~,~,l]=unique({d(:).cellType});

q = false(numel(allT),1);


for f=1:numel(folds)
    fid = fopen(sprintf(fname,f),'wb');
    
    %save number of cells as int
    fwrite(fid,numel(folds{f}),'uint16');
    fwrite(fid,numel(allT),'uint16');
    fwrite(fid,numel(allS),'uint16');
    
    
    %start by saving the start time, end time, and dt as double
    fwrite(fid,allT([1 end]),'double');
    fwrite(fid,dt,'double');
    
    %save smallest and largest spot sizes as int
    fwrite(fid,allS([1 end]),'uint16');
    
    %for each cell
    d2 = d(folds{f});
    l2 = l(folds{f});
    minT2 = minT(folds{f});
    maxT2 = maxT(folds{f});
    se2 = se(folds{f});
    allD2 = allD(folds{f},:,:);
    for n=1:numel(folds{f})
        
        %save length of cellName as int, then cellName as char
        if isa(d2(n).cellName,'cell')
            fwrite(fid,length(d2(n).cellName{1}),'uint8');
            fwrite(fid,d2(n).cellName{1},'char');
        else
            fwrite(fid,length(d2(n).cellName),'uint8');
            fwrite(fid,d2(n).cellName,'char');
        end
        
        %save cellType as char (1 to 32)
        fwrite(fid,l2(n),'char');
        
        %save quality for each time point as logical
        q(:) = false;
        q(allT>=minT2(n) & allT<=maxT2(n)) = true;
        fwrite(fid,q,'logical');
        
        %save quality for each spot size...
        fwrite(fid,uint16(se2{n}),'uint32');
        %save rebinned data
        fwrite(fid,squeeze(allD2(n,:,:))','double');
        
    end
    
    fclose(fid);
end