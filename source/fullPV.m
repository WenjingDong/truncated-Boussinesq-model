function [q0, q1, q2, q0linear, q1linear, q2linear] = fullPV(u0, v0, u1, v1, b1, u2, v2, b2, N, f, m, Lx)

% set up wavenumber
nx = size(u0,1); ny = size(u0,2);


% test
%u0 = 0 * u0;
%v0 = 0 * v0;
%u2 = 0 * u2;
%v2 = 0 * v2;
%b2 = 0 * b2;


n = nx/2;  
kmax = n - 1;     
[kx_,ky_] = ndgrid([0:kmax -kmax-1:-1],[0:kmax -kmax-1:-1]);
kx_ = kx_ * 2*pi/Lx;
ky_ = ky_ * 2*pi/Lx;

ikx_ = sqrt(-1) * kx_;
iky_ = sqrt(-1) * ky_;



    u0k = fft2(u0)/(nx*ny);
    v0k = fft2(v0)/(nx*ny);
    u1k = fft2(u1)/(nx*ny);
    v1k = fft2(v1)/(nx*ny);
    w1k = 0 * -1/m * (ikx_ .* u1k + iky_ .* v1k); % FOR HYDROSTATIC PV
    
    b1k = fft2(b1)/(nx*ny);
    u2k = fft2(u2)/(nx*ny);
    v2k = fft2(v2)/(nx*ny);
    b2k = fft2(b2)/(nx*ny);
    w2k = 0 * -1/(2*m) * (ikx_ .* u2k + iky_ .* v2k); % FOR HYDROSTATIC PV
    
    v0x = ifft2(ikx_ .* v0k, 'symmetric') * (nx*ny);
    v1x = ifft2(ikx_ .* v1k, 'symmetric') * (nx*ny);
    b1x = ifft2(ikx_ .* b1k, 'symmetric') * (nx*ny);
    v2x = ifft2(ikx_ .* v2k, 'symmetric') * (nx*ny);
    b2x = ifft2(ikx_ .* b2k, 'symmetric') * (nx*ny);
    w1x = ifft2(ikx_ .* w1k, 'symmetric') * (nx*ny);
    w2x = ifft2(ikx_ .* w2k, 'symmetric') * (nx*ny);
    
        
    u0y = ifft2(iky_ .* u0k, 'symmetric') * (nx*ny);
    u1y = ifft2(iky_ .* u1k, 'symmetric') * (nx*ny);
    b1y = ifft2(iky_ .* b1k, 'symmetric') * (nx*ny);
    u2y = ifft2(iky_ .* u2k, 'symmetric') * (nx*ny);
    b2y = ifft2(iky_ .* b2k, 'symmetric') * (nx*ny);
    w1y = ifft2(iky_ .* w1k, 'symmetric') * (nx*ny);
    w2y = ifft2(iky_ .* w2k, 'symmetric') * (nx*ny);
    
    % full PV
    q0 = N^2 * (v0x - u0y) + ...
    1/2 * ( b1x .* (w1y +   m * v1) - b1y .* (  m * u1 + w1x) +   m * b1 .* (v1x - u1y) + ...
         0 *  b2x .* (w2y + 2*m * v2) - 0* b2y .* (2*m * u2 + w2x) +  0* 2*m * b2 .* (v2x - u2y));
        
    q1 = f * m * b1 + N^2 * (v1x - u1y) + m * b1 .* (v0x - u0y) ...
     +   1/2  * ( b2x .* (w1y + m * v1) +    b1x .* (w2y + 2 * m * v2)...
                     -b2y .* (w1x + m * u1) -    b1y .* (w2x + 2 * m * u2)...
            +2 * m * b2  .* (v1x - u1y)    + m * b1 .* (v2x - u2y) );
          
    q2 =  f * 2 * m * b2 +  N^2 * (v2x - u2y) + 0* 2 *m *b2 .* (v0x - u0y)...
    +   1/2 * ( -b1x .* (w1y + m * v1) +  b1y .* (w1x + m * u1) + m * b1 .* (v1x - u1y) );
     
    % PV, neglecting mode 0, mode 2 interaction and mode 1, mode 2, mode 0
    % interaction
    
    %q0 =  N^2 * (v0x - u0y) + ...
    %1/2 * (  b1x.*( w1y+m*v1 ) -  b1y.*( m*u1+w1x ) + m*b1.*( v1x - u1y ) );
 
    %q1 = f * m * b1 + N^2 * (v1x - u1y);
    
    %q2 = f * 2*m * b2 + N^2 * (v2x - u2y) ...
    %+ 1/2*( -b1x.*(w1y + m*v1) + b1y.*(w1x + m*u1) + m*b1.* (v1x - u1y) );
    
    q0linear =  N^2 * (v0x - u0y);
    q1linear = f * m * b1 + N^2 * (v1x - u1y);
    q2linear = f* 2 * m * b2 + N^2 * (v2x - u2y);



end
