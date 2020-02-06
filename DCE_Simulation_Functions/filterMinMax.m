function filterMinMax(inFileName,outFileName,lowerBound,upperBound)
%Filter image values so that only nans and values within bounds are
%preserved. Values outside range are replaced by lowerBound or upperBound.
%inFileName: input path and filename
%outFileName: output path and filename for filtered image
%lowerBound
%upperBound

[path,name,ext]=fileparts(inFileName);

inFileHdr=spm_vol(inFileName);
[inSI,temp]=spm_read_vols(inFileHdr);

outSI=inSI;
outSI(inSI<=lowerBound)=lowerBound;
outSI(outSI>=upperBound)=upperBound;

outFileHdr=inFileHdr;
outFileHdr.fname=outFileName;
spm_write_vol(outFileHdr,outSI);

end