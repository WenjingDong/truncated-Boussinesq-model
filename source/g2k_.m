function [fk] = g2k_(fg,dealias)

nx = size(fg,1);
ny = size(fg,2);


normfac = 1/(nx*ny);
fbk = normfac*fft(fft(fg,[],1),[],2);
if dealias, 
    n = nx/3;
    fk = [fbk(1:n,1:n,:); fbk(2*n+1:end,1:n,:)];
else
    n = nx/2;
    fk = fbk(:,1:n,:);
end

fk(n+1:end,1,:) = 0;  % this part given by C-S
%fk(1,1,:) = 0;
end