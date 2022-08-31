function [fg] = k2g_(fk,dealias)
    
normfac = ( size(fk,1) )^2;
if dealias, normfac = (3/2)^2*normfac; end
fg = normfac*real(ifft(ifft(sympad_(fk,dealias),[],1),[],2));
end

