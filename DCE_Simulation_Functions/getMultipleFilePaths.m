function [paths,names]=getMultipleFilePaths(str)
%used to get full paths and names of files specified using a wildcard
%returns full paths of file
%returns error if there are no matching files

targetFiles=dir(str);
if size(targetFiles,1)<1; disp('No matching files found.'); paths={}; names={}; return; end;

targetDir=str(1:max(strfind(str,'/')));

paths=cell(size(targetFiles,1),1);
names=cell(size(targetFiles,1),1);

for iFile=1:size(targetFiles,1)
    paths{iFile}=[targetDir targetFiles(iFile).name];
    names{iFile}=targetFiles(iFile).name;
end

end