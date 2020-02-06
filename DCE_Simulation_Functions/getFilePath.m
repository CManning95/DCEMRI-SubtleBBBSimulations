function [path,name]=getFilePath(str)
%used to get full path of a single file whose name is not known but can be
%specified using a wildcard
%returns full path of file
%returns error if there is more than one possible file

targetFile=dir(str);
if size(targetFile,1)~=1; error('There is not exactly 1 target file specified.'); end;
targetDir=str(1:max(strfind(str,'/')));
path=[targetDir targetFile(1).name];
name=targetFile(1).name;

end