function ROIData=measure4D(image,ROI)
%ROIData: struct containing mean, median etc. for ROI; each variable is a
%column vector containing value at each time point
%image, mask: either NIFTI files containing 4D image and ROI (if variable is a string), or the corresponding
%image matrices

if isstr(image)
    [signal,temp]=spm_read_vols(spm_vol(image)); %read file containing values of interest
else
    signal=image;
end

if isstr(ROI)
    [mask,temp]=spm_read_vols(spm_vol(ROI)); %read mask file
else
    mask=ROI;
end

if ~isempty(find(mask ~= 1 & mask ~= 0)); error(['Mask contains voxels which are not 0 or 1!']); end; %check masks contain only 0 and 1

N=size(signal,4);

ROIData.mean=nan(N,1);
ROIData.median=nan(N,1);
ROIData.SD=nan(N,1);

for n=1:N %loop through time points of 4D image
    signalThisTimepoint=squeeze(signal(:,:,:,n));
    ROIData.mean(n,1) = nanmean(signalThisTimepoint(mask==1));
    ROIData.median(n,1) = nanmedian(signalThisTimepoint(mask==1));
    ROIData.SD(n,1) = nanstd(signalThisTimepoint(mask==1));
end

end
