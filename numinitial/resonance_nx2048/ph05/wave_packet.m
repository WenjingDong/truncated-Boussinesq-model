% wave packet test, non-dimensionalized by 1/k
clear;

load('input.mat');

% dimensional physical parameters
H_dim = 4000;
m_dim = pi/H_dim;
N_dim = 2e-3;
f_dim = 1e-4;
k_dim = 7* f_dim/N_dim  * m_dim; % comovig solution if less than sqrt(1/3) *f_dim/N_dim * m_dim



% non-dimensional physical parameters
k = 1;
m = m_dim/k_dim * k;

N = 1;
f = f_dim/N_dim * N;

omega = sqrt(N^2 * k^2/m^2 + f^2);
cg = N^2/m^2 * k/omega;
cp = omega/k;
c = N/(2*m);

mu = 0.05;

if (c - cg<0)
    disp('no comoving solution');
end

% initial conditions, non-dimensional
alpha = 8 * round(1/mu); % number of wave lengths in the domain.
nx =  2048; 
ny = nx;
Lx = alpha * 2*pi/k;
x = linspace(0,Lx * (nx-1)/nx,nx) - Lx/16;
y = linspace(0,Lx * (ny-1)/ny,ny) - Lx/2;
[x_,y_] = ndgrid(x,y);
dx = Lx/nx;

% wavenumber matrix
[kx_,ky_] = ndgrid([0:nx/2-1 -nx/2:1:-1], [0:ny/2-1 -ny/2:1:-1]); 
kx_ = 2*pi/Lx*kx_; ky_ = 2*pi/Lx * ky_;
ikx_ = sqrt(-1)*kx_; 
iky_ = sqrt(-1)*ky_;

omega_ = sqrt(N^2* (kx_.^2 + ky_.^2)/m^2 + f^2);
omega_(nx/2+1:end,:) = -omega_(nx/2+1:end,:);
omega2_f2 = omega_.^2 - f^2;
omega2_f2(1,1) = 0.1; % avoid divide zero

AMP = 0.01;
a = AMP *exp(-mu^2*k^2*(x_.^2 + y_.^2)/2);
da2dx = AMP^2 * (-2*mu^2*k^2*x_).* exp(-mu^2*k^2*(x_.^2 + y_.^2));
da2dy = AMP^2 * (-2*mu^2*k^2*y_).* exp(-mu^2*k^2*(x_.^2 + y_.^2));

% vertical velocity for phase 1
theta = pi/2;
w = omega/m * a .* cos(k*x_ + theta);

% hydrostatic
what = fft2(w)/(nx*ny);
what(1,1) = 0;
uhat = 1./omega2_f2.* N^2./omega_.^2.*omega_./m .* (ikx_.*omega_ - ky_*f).*what;
vhat = 1./omega2_f2.* N^2./omega_.^2.*omega_./m .* (iky_.*omega_ + kx_*f).*what;
bhat = -sqrt(-1) * N^2./omega_ .* what;

u1 = nx * ny * ifft2(uhat,'symmetric');
v1 = nx * ny * ifft2(vhat,'symmetric');
b1 = nx * ny * ifft2(bhat,'symmetric');

b2 = N^2/(4*m)*a.^2;

Sin = zeros(nx,nx,8);

Sin(:,:,3) =  u1/sqrt(2);
Sin(:,:,4) =  v1/sqrt(2);
Sin(:,:,5) =  b1/sqrt(2);
%
Sin(:,:,1) = u0;
Sin(:,:,2) = v0;
Sin(:,:,8) = b2/sqrt(2);

vel_ind = [1 2 3 4 6 7];
umax = max([max(max(max(abs(Sin(:,:,vel_ind))))) cp N/(2*m)]);



% time step
cfl = 0.05;
dt = cfl * dx/umax;
 
T =   7/8*Lx/cg;
numbersteps = floor(T/dt);
savesteps = floor(numbersteps/20);
writeoutput = 1;
tic
[Sout,Diag,Diag2] = bous3m_wd2(Sin, N, f, m, numbersteps, savesteps, dt, Lx, writeoutput);
 toc;
 
time = Diag2.time;
ke0 = Diag2.ke0; ke1 = Diag2.ke1; ke2 = Diag2.ke2;
pe1 = Diag2.pe1; pe2 = Diag2.pe2;



save('output.mat','Diag','Diag2');
