function write_field(field,fname,frame,ftype,dtype,unit)

% WRITE_FIELD(field,fnm,frame,ftype,dtype,unit)  
%     Writes 'field' 
%     to 'frame' of direct access binary file with name 'fname'.
%     Optional 'ftype' and 'dtype' have defaults 's' and 'real*8', 
%     respectively.  Optional 'unit' is byte-width of each element
%     of type dtype (default is 8 which is correct for real*8).
%
%     See also READ_FIELD, FOPEN, FWRITE.

% Defaults:
unitd = 8; dtyped = 'real*8'; ftyped = 'n'; framed = 1;

switch nargin
   case 2, unit = unitd; dtype = dtyped; ftype = ftyped; frame = framed;
   case 3, unit = unitd; dtype = dtyped; ftype = ftyped;
   case 4, unit = unitd; dtype = dtyped; 
   case 5, unit = unitd;
end

nx = size(field,1); ny = size(field,2); nz = size(field,3);

fnm = strcat(fname,'.bin');
[fid,msg]=fopen(fnm,'r',ftype);
if fid>=0
  disp('WARNING: output file exists -- field will be appended')
else
  disp('Output file does not exist -- will create it')
end
[fid,msg]=fopen(fnm,'a',ftype);
if fid<0, disp(fnm), error(msg), end

switch isreal(field)
case 1
   status = fseek(fid,unit*nx*ny*nz*(frame-1),-1);
   if status~=0 disp(ferror(fid)), fclose(fid), return, end  
   fwrite(fid,field,dtype);
   fclose(fid);
case 0         % its a spectral field with an imaginary part
   status = fseek(fid,unit*nx*ny*nz*2*(frame-1),-1);
   if status~=0 disp(ferror(fid)), fclose(fid), return, end  
   disp(status)
   status = fwrite(fid,real(field),dtype);
   disp(status)
   status = fwrite(fid,imag(field),dtype);
   disp(status)
   fclose(fid);
end

