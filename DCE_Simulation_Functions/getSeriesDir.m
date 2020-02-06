function seriesDirStr=getSeriesDir(parentDirStr,seriesNum)
%find directory beginning with seriesNum or seriesNum_
%useful for getting directories containing a scan series

temp=dir([parentDirStr '/']);
temp2=~cellfun(@isempty,regexp({temp.name},...
    ['^' num2str(seriesNum) '_']))... %look for directories starting with seriesNum_
    | strcmp({temp.name},num2str(seriesNum)); %look for directories with names seriesNum

if sum(temp2)~=1; error('Cannot find single unique dicom directory for this series.'); end

seriesDirStr=[parentDirStr '/' temp(temp2).name];

end