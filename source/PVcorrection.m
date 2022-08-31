% This function computes u0, v0, and b2 to make PV uniform
function [u0, v0, b2] = PVcorrection(uhat, vhat, what, bhat, u1, v1,  b1, ikx_, iky_, N, f,  m)
   [nx, ny] = size(ikx_);
   
   IK2 = ikx_.^2 + iky_.^2;
   IK2(1,1) = 0.1;
   
   v1x = nx * ny * ifft2(ikx_ .* vhat, 'symmetric');
   b1x = nx * ny * ifft2(ikx_ .* bhat, 'symmetric');
   w1x = nx * ny * ifft2(ikx_ .* what, 'symmetric');

        
   u1y = nx * ny * ifft2(iky_ .* uhat, 'symmetric');
   b1y = nx * ny * ifft2(iky_ .* bhat, 'symmetric');
   w1y = nx * ny * ifft2(iky_ .* what, 'symmetric');
   
   b2 = -1/2 * ( -b1x .* (w1y + m * v1) + b1y .* (w1x + m * u1) + m * b1 .* (v1x - u1y) );
   b2 = b2/(f * 2 * m);
 
   
   rhspsi0 = -1/2 * ( b1x .* (w1y + m * v1) - b1y .* ( m * u1 + w1x) + m * b1 .* (v1x - u1y) );
   rhspsi0hat = fft2(rhspsi0)/(nx * ny);
   rhspsi0hat(1,1) = 0; % force it to be zero,  corresponds to a constant
   psi0hat = 1/N^2 ./IK2 .* rhspsi0hat;
   
   
 
   
   v0hat =  ikx_ .* psi0hat;
   u0hat = -iky_ .* psi0hat;
  
   u0 = nx * ny *  ifft2(u0hat,'symmetric');
   v0 = nx * ny *  ifft2(v0hat,'symmetric');
   
   
 
end
        