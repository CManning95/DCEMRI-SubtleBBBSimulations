function SPMWrite4D(vol4D,data,outDir,fileName,datatype)
%write 3D or 4D file using SPM header struct array as template
%vol4D: template 3D or 4D SPM header struct array
%data: 3D or 4D data
%outdir: output directory
%filename: without extension
%datatype: spm datatype, e.g. 16 for float

vol3D=vol4D(1);
NFrames=size(data,4);

for iFrame=1:NFrames
    volOut=vol3D;
    volOut.dt=[datatype 0];
    volOut.fname=[outDir '/' fileName '_temp_' num2str(iFrame,'%02d') '.nii'];
    spm_write_vol(volOut,data(:,:,:,iFrame));
end

spm_file_merge(sort(getMultipleFilePaths([outDir '/' fileName '_temp_*.nii'])),[outDir '/' fileName '.nii'],0);

delete([outDir '/' fileName '_temp_*.nii']);

end
