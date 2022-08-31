function [hkb] = sympad_(hk,dealias)

% Input hk is assumed to be upper-half-plane spectral field, size
% [nx,nx/2], representing Fourier transform of real field.  Use
% conjugate-symmetry to make hbk with size [nx,nx].  Optional
% argument 'dealias', if present and true, creates [3*nx/2,3*nx/2]
% padded field for making dealiased product
    nx = size(hk,1);

    n = nx/2;
    if nx~=2*n, error('need size(hk,1) = 2*size(hk,2)'), end 
    
    p = 0;
    if dealias, p = 1; end

    hkb = zeros((2+p)*n,(2+p)*n,size(hk,3));        % Full-size array
    
    hk(n+2:end,1,:) = conj(hk(n:-1:2,1,:));         % Make ky=0 C-S
    hkb(1:n,1:n,:) = hk(1:n,:,:);                   % [0..K,:]  
    hkb((1+p)*n+2:end,1:n,:) = hk(n+2:end,:,:);     % [-K..-1,:]
    
    hkb(2:end,(1+p)*n+2:end,:) = conj(hkb(end:-1:2,n:-1:2,:)); 
    hkb(1,(1+p)*n+2:end,:) = conj(hkb(1,n:-1:2,:));
 
end
    
   