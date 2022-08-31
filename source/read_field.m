function field = read_field(file,nx,ny,nz,frmvec,is_real,ftype,dtype,unit)

% field = READ_FIELD(file,nx,ny,nz,frmvec,is_real,ftype,dtype,unit)
%     Read direct access 'file'//.bin.  All parameters except 'file'
%     are optional.  If specified, [nx ny nz] is assumed as 
%     size of 1 frame of 'field'.  If no options are specified, 
%     file is assumed to be a 0-d series, and the entire file 
%     is read into first dimension of 'field'.  If nx == 2*ny-1, 
%     field is assumed to be spectral, and 'file'//.bin
%     is expected to contain staggered frames of real and imaginary 
%     components  of field.  If 'frmvec' is specified, frames in this 
%     vector are read.  In the case of 1-d fields [nx 1 1], 
%     frames are stored in 2nd dimension.  In 2-d case, frames are 
%     stored in 3rd dimension.  In 3-d case, they are stored 
%     successively in 3rd dimension.
%  
%     Optional 'ftype', 'dtype', and 'unit' specify, respectively, 
%     binary storage protocol (e.g. IEEE-Big Endian), data type (e.g.
%     real*8) and the number of bytes in a record for the file.
%     Current defaults are ftype = 'ieee_be.l64', 
%     dtype = 'real*8', and unit = 8.
%     Optional 'is_real' can be used (set to 1) to override default assumption
%     that input field is imaginary (which is assumed if nx =
%     2*ny-1, because that is the relationship for spectral outputs
%     from SQG model).
%
%     See also FOPEN, FREAD, WRITE_FIELD.

% Defaults:
unitd = 8; dtyped = 'real*8'; 
% ftyped = 's';   % to read in Max files
ftyped = 'n';    
framed = 1; nzd = 1; nyd = 1; nxd = 1;

% Assume its a complex field if dimensions are those of complex
% field from SQG
if (nargin>2 & nx==2*ny-1)
  is_reald=0;
else
  is_reald=1;
end

switch nargin
  case 1, unit = unitd; dtype = dtyped; ftype = ftyped; is_real = is_reald;
          frmvec = framed; nz = nzd; ny = nyd; nx = nxd;
  case 2, unit = unitd; dtype = dtyped; ftype = ftyped; is_real = is_reald; 
          frmvec = framed; nz = nzd; ny = nyd;
  case 3, unit = unitd; dtype = dtyped; ftype = ftyped; is_real = is_reald; 
          frmvec = framed; nz = nzd;
  case 4, unit = unitd; dtype = dtyped; ftype = ftyped; is_real = is_reald;
          frmvec = framed; 
  case 5, unit = unitd; dtype = dtyped; ftype = ftyped; is_real = is_reald;
  case 6, unit = unitd; dtype = dtyped; ftype = ftyped; 
  case 7, unit = unitd; dtype = dtyped; 
  case 8, unit = unitd;
end

% Open file
fnm = strcat(file,'.bin'); 
[fid,msg] = fopen(fnm,'r',ftype);
if (fid<0) 
  disp(fnm)
  disp(msg)
  field = 0;
  return
end

switch nx
case 1     % Time series or rank-1 array (vector)
   field = fread(fid,[nx inf],dtype);
otherwise   
   j=0;
   for frm = frmvec
      j=j+1;
      yvec = (j-1)*ny*nz+1:j*ny*nz;
      if (is_real==0)
         fseek(fid,unit*nx*ny*nz*2*(frm-1),-1);
         fieldr = fread(fid,[nx ny*nz],dtype);
         fieldi = fread(fid,[nx ny*nz],dtype);
         field(1:nx,yvec) = fieldr + i*fieldi;
      else
         fseek(fid,unit*nx*ny*nz*(frm-1),-1);
         field(1:nx,yvec) = fread(fid,[nx ny*nz],dtype);
      end
    end
    
    % Layerize
    temp = field;
    field = zeros(nx,ny,nz,length(frmvec));
    for frm = 1:length(frmvec)
      for layer = 1:nz
	lindi = ((frm-1)*nz+layer-1)*ny+1;
	uindi = ((frm-1)*nz+layer)*ny;
	field(:,:,layer,frm) = temp(:,lindi:uindi);
      end
    end
    field = squeeze(field);  % in case nz=1 or frames=1;    

end

fclose(fid);   % Close file
