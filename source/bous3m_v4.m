function [Sout,Diag,Diag2,q] = bous3m_v4(Sin,N,f,m, numsteps,savestep,dt,Lx,writeoutput)

%  [Sout,Diag,hmov] = bous3m(S,Bu,Ro,m,numsteps,savestep)
%
%  Solves truncated HS Boussinesq model:  BT + BCm + BC2m 
%  Domain:  (x,y) \in [0,2pi) and z in [-pi,0]
%
%  Inputs
%
%  Sin:       N x N x 7 array with initial zeta0, u1, v1, b1,
%                 u2, v2, b2, respectively
%  Ro:        U/(f*L), where dimensional (x,y) \in [0,2*pi*L]
%  Bu:        (N*H/(f*L))^2, where dimensional z \in [-H*pi,0]    
%  m:         Vertical wavenumber of BC1 wave (then BC2 wave has
%                 2*m). Must be integer. 
%  numsteps:  Total number of timesteps 
%  savestep:  Frequency, in timesteps, to save output    
%
%  Outputs
%
%  Sout:      Arranged as S, but with 4th dimension for time
%  Diag:      Struct containing diagnostics 
%    
%  Numerical details
%
%  Model is spectal, in square domain of size 2*pi x 2*pi X pi.  Input
%  fields must have nx = ny = 2^(integer).  Nonlinear terms are done in
%  physical space using dealiased product.  Uses AB3 timestepping with
%  trapezoidal hyperviscosity of order a.  Timestep and hyperviscosity
%  are set adaptively, via dt = dttune*dx/max(|u|) and nu =
%  nutune*dx/(dt*kmax^a), where dx = 2*pi/nx, and kmax = nx/2 - 1.
%
%  Tuning factors dttune and nutune and hyperviscous order a
%  can be set by editing this file directly.
%  
%  TO DO:  use varargin param, value list or input param struct
%          switch transforms to natural FFT k-grid
%          streamline RHS function with ufunc on structs
%          test removing loops over dim 3 in fft2 and ifft2
    
hvord  = 8;         % hyperviscosity order, nu del^hvord u 
nutune = 0.5;        % tuning for adaptive hyperviscosity orginally 0.5

    
% Get and check dimensions
[nx,ny,nfields] = size(Sin);
    
%if (nx~=ny), error('must have nx = ny'), end
%if (mod(log2(nx),1)~=0), error('must have nx = 2^n, n integer'), end
%if (nfields~=8), error('input state S must have size [nx,nx,7]'), end

% Check for outputs requested

%if (nargout>2), makemov=true;  end

% Named indeces to select fields from S, e.g. S(:,:,u1_ind) is u1(:,:)
u0_ind = 1;  v0_ind = 2;
u1_ind = 3;  v1_ind = 4;  b1_ind = 5; 
u2_ind = 6;  v2_ind = 7;  b2_ind = 8;

% Remove the following fields from calculation (for testing)
kill_list = [];

% Only filter (hyperdiffuse) the following fields
filt_field_list = [u0_ind v0_ind u1_ind v1_ind u2_ind v2_ind]; 
%filt_field_list = [];
% Domain width = 2*pi, so:
dx = Lx/nx;             

% For RHS 
sqr2 = sqrt(2); 
sqr2i = 1/sqrt(2);

% For use in calls to k2g and g2k
da = true;           
noda = false;         

% Max velocity for CFL check 
%Cmax = sqrt(Bu/m^2 + 1);  % Max phase speed (omega^2 = Bu*K^2/m^2+1)
Cmax = 0;

% Set up arrays of wavenumbers for computing derivatives.  Spectral
% arrays (ending with 'k') are size [2*n,n], where n = nx/2, which
% omits redundant -ky part

n = nx/2;  
kmax = n - 1;     
[kx_,ky_] = ndgrid([0:kmax -kmax-1:-1],0:kmax);
kx_ = kx_ * 2*pi/Lx;
ky_ = ky_ * 2*pi/Lx;

K_   = sqrt(kx_.^2 + ky_.^2);
ikx_ = 1i*kx_; 
iky_ = 1i*ky_;
IK2_ = K_.^(-2);

% This is for computation of BT velocities u0, v0 from BT vorticity
% z0.  In the BT vorticity forumlation, we've removed the
% homogeneous mode, and want to prevent divide-by-zero, so we set
% the k,l = 0,0 element of 1/K^2 to 0. 
% NOTE:  do we need to keep track of that neglected homogeneous
% part?  
IK2_(1,1) = 0;  

% Initialize fftw
fftw('dwisdom',[]);
fftw('planner','measure');

% Get initial spectral fields
Sk = g2k(Sin,noda); % nx * n * nfields

% Use best fftw 'plan' after first use of fft/ifft (google matlab fftw)
fftinfo = fftw('dwisdom');
fftw('dwisdom',fftinfo);

    % Trapezoidal hyperdiffusion operators 
nudto2 = nutune*dx/(2*kmax^(hvord-1));         % nudto2 = nu*dt/2
nu = 2*nudto2/dt;
    
filtL = ones(size(Sk));
filtR = ones(size(Sk));
for j = filt_field_list
    filtR(:,:,j) = (1+nudto2*K_.^hvord).^(-1);
    filtL(:,:,j) = (1-nudto2*K_.^hvord).*filtR(:,:,j);
end   

clear K_ % not needed after above.

% declare variables: 3 RHS stage arrays and output fileds
Rk = zeros(size(Sk,1),size(Sk,2),size(Sk,3),3);

% Preallocate output arrays
nframes = floor(numsteps/savestep)+1;   % number of output fields
Diag.time = zeros(1,nframes);       % model time of output 
Diag.ke = zeros(1,nframes);         % kinetic energy
Diag.pe = zeros(1,nframes);         % potential energy

if writeoutput
   Sout = zeros(nx,nx,nfields);
   q = zeros(nx,nx,3);
else
   Sout = zeros(nx,nx,nfields,nframes);
   q = zeros(nx,nx,3,nframes);
end

% Initialize output arrays with starting fields
%kss:  These diagnostics are all computed at first timestep in
%      (mod(counter-1,savestep)==0) in main loop, right?

%Sout(:,:,:,1) = k2g(Sk,noda);
%Diag.time(1) = 0;
%Diag.ke(1) = mean(mean(Sout(:,:,u0_ind,1).^2)) + mean(mean(Sout(:,:,v0_ind,1).^2)) + ...
%    mean(mean(Sout(:,:,u1_ind,1).^2)) + mean(mean(Sout(:,:,v1_ind,1).^2)) + ...
%    mean(mean(Sout(:,:,u2_ind,1).^2)) + mean(mean(Sout(:,:,v2_ind,1).^2));
%Diag.ke(1) = Diag.ke(1)/2;
%Diag.pe(1) = 1/N^2 * (mean(mean(Sout(:,:,b1_ind,1).^2)) + mean(mean(Sout(:,:,b2_ind,1).^2)))/2;

Diag2.time = zeros(numsteps,1);
Diag2.ke = zeros(numsteps,1);
Diag2.pe = zeros(numsteps,1);
Diag2.ke0 = zeros(numsteps,1);
Diag2.ke1 = zeros(numsteps,1);
Diag2.ke2 = zeros(numsteps,1);
Diag2.keu2 = zeros(numsteps,1);
Diag2.kev2 = zeros(numsteps,1);
Diag2.pe1 = zeros(numsteps,1);
Diag2.pe2 = zeros(numsteps,1);

x = linspace(0,Lx*(nx-1)/(nx),nx) - Lx/4;
y = linspace(0,Lx*(ny-1)/(ny),ny) - Lx/2;
[x_,y_] = ndgrid(x,y);      

% Set counters
counter = 1;  % time step index
frame = 1;    % diagnostic output counter % 1st frame is initial value
t = 0;        % model time

keepgoing = true;

while keepgoing
       
    % Check energy conservation by first cutting out modes
    Sk(:,:,kill_list) = 0;
    
    Umax = getrhs(counter, frame);  % Sets Rk and q
     
    % Exit if blowing up, and save last field.
    if (Umax>1e6 | isnan(Umax))  
        disp(strcat('Blow up! Umax= ',num2str(Umax),', t= ',num2str(t)))
        Diag.time(frame+1) = t;
	save
        keepgoing = false;
    end
    
     % Save output at frequency savestep
    if (mod(counter-1,savestep)==0)  
     
        if writeoutput 
            %save(outputfile,strcat('Sout',num2str(frame)),'-append')
 	    Sout = k2g(Sk,noda);  
	    write_field(Sout,'Sout',frame);
	   % write_field(q,'PV',frame);
	    write_field(t,'time',frame);
	    %save
	    frm = 1;
	else
	    Sout(:,:,:,frame) = k2g(Sk,noda);  
	    qout(:,:,:,frame) = q;
	    frm = frame;
        end

        Diag.time(frame) = t;
        
        div = sum(sum(k2g(ikx_.* Sk(:,:,u0_ind) + iky_.*Sk(:,:,v0_ind), noda)));
        KE0  = mean(mean(Sout(:,:,u0_ind,frm).^2 + Sout(:,:,v0_ind,frm).^2))/2;
        KE12 = mean(mean(Sout(:,:,u2_ind,frm).^2 + Sout(:,:,v2_ind,frm).^2 +...
                         Sout(:,:,u1_ind,frm).^2 + Sout(:,:,v1_ind,frm).^2))/2;
        PE = 1/N^2 * mean(mean(Sout(:,:,b1_ind,frm).^2 + Sout(:,:,b2_ind,frm).^2))/2;
        Diag.ke(frame) = KE0  + KE12;
        Diag.pe(frame) = PE;
        
        %disp(strcat('sum(psi J) =',num2str(sum(sum(psi.*foo)))))
        disp(strcat('Wrote frame >',num2str(frame),' out of >',num2str(nframes)))
        disp(strcat('max(|u|) = ',num2str(Umax),', dt = ',num2str(dt),', nu = ',num2str(nu)))
        disp(strcat('KE1 + KE2 = ',num2str(KE12),', PE = ', num2str(PE), ',sum E = ', num2str(KE0+KE12+PE)));
        disp(strcat('div u0 = ', num2str(div)));        

        % update frame
        frame = frame + 1;
    end
    
    % Timestep
    Sk = filtL.*Sk;
    if (counter==1)  % Euler step
        Sk = Sk + filtR.*(dt*Rk(:,:,:,1));
    end
    if (counter==2) % AB2 step
        Sk = Sk + filtR.*((dt*3/2)*Rk(:,:,:,1) - (dt/2)*Rk(:,:,:,2));
    end
    if (counter>2) % AB3
        Sk = Sk + filtR.*((dt*23/12)*Rk(:,:,:,1) - (dt*16/12)*Rk(:,:,:,2) + (dt*5/12)*Rk(:,:,:,3));
    end 
       
    % update counter and time
    counter = counter + 1;
    t = t+dt;  % clock
   
    if (counter-1==numsteps), disp('End reached'), keepgoing=false; end
    
end % timestepping loop

%-------------------------------------------------------------------
% Internal functions: function end at end of file means all
%                     internal functions have access to all
%                     variables in main  
%-------------------------------------------------------------------


function [Umax] = getrhs(counter, frame)
   
    persistent S p0k w1 w2 w1k w2k Lk NLx NLy NLz NLxy 
    % persistent effectively allocates fixed memory, upon first
    % use, for temp variables accessed only in this function
    % NOTE:  all gridded fields here are size 3/2*[nx,nx]

    % All 7 equations of form:  d_t u = L - Ro(d_x NLx + d_y NLy + m NLz)
    % All pressure terms from buoyancy:  p_j = -b_j/j, j = m or 2*m

    % Save past stages for AB3
    Rk(:,:,:,3) = Rk(:,:,:,2);
    Rk(:,:,:,2) = Rk(:,:,:,1);

    % Get gridded fields for NL products
    S = k2g(Sk,da);  
    
    % w_j = -(d_x u_j + d_y v_j)/j, where j = m or 2*m; 
    
    w1k = -(ikx_.*Sk(:,:,u1_ind)+iky_.*Sk(:,:,v1_ind))/m;
    w2k = -(ikx_.*Sk(:,:,u2_ind)+iky_.*Sk(:,:,v2_ind))/(2*m);
    w1 = k2g(w1k,da);
    w2 = k2g(w2k,da);
    
    % ----- RHS p0 ------
    NLx =g2k(S(:,:,u0_ind).^2 + S(:,:,u1_ind).^2 + S(:,:,u2_ind).^2, da);

    NLxy = g2k(S(:,:,v0_ind).*S(:,:,u0_ind) + S(:,:,v1_ind).*S(:,:,u1_ind) + ...
                S(:,:,v2_ind).*S(:,:,u2_ind), da);
    
    NLy = g2k(S(:,:,v0_ind).^2 + S(:,:,v1_ind).^2 + S(:,:,v2_ind).^2, da);
    
    Lk = f* (iky_.* Sk(:,:,u0_ind) - ikx_.* Sk(:,:,v0_ind));
    
    p0k = Lk + (ikx_.^2 .* NLx + iky_.^2.* NLy + 2*ikx_ .*iky_ .* NLxy);
    p0k = p0k.*IK2_;
    
    % ----- RHS u0, v0 ------
    Rk(:,:,u0_ind,1) =  f*Sk(:,:,v0_ind) - ikx_.* p0k -  ikx_.* NLx  - iky_.*NLxy ;
    Rk(:,:,v0_ind,1) = -f*Sk(:,:,u0_ind) - iky_.* p0k -  ikx_.* NLxy - iky_.*NLy;
    
    % ----- RHS u1 ----- 
    NLx = g2k(2*S(:,:,u1_ind).*S(:,:,u0_ind) + sqr2*S(:,:,u1_ind).*S(:,:,u2_ind) ,da);
    NLy = g2k( S(:,:,v1_ind).*S(:,:,u0_ind) + S(:,:,u1_ind).*S(:,:,v0_ind) ...
               + sqr2i*(S(:,:,u1_ind).*S(:,:,v2_ind) + S(:,:,u2_ind).*S(:,:,v1_ind)) ,da);
    NLz = g2k( S(:,:,u0_ind).*w1 + sqr2i*(S(:,:,u1_ind).*w2 - S(:,:,u2_ind).*w1) ,da);

    Lk  = f * Sk(:,:,v1_ind) + 1/m*ikx_.*Sk(:,:,b1_ind);

    Rk(:,:,u1_ind,1) = Lk - (ikx_.*NLx + iky_.*NLy + m*NLz);
     
    %rhsu1 =  k2g(Rk(:,:,u1_ind,1),noda);

    % ----- RHS v1 ----- 
    NLx = g2k(S(:,:,u1_ind).*S(:,:,v0_ind) + S(:,:,u0_ind).*S(:,:,v1_ind) ...
              + sqr2i*(S(:,:,u2_ind).*S(:,:,v1_ind) + S(:,:,u1_ind).*S(:,:,v2_ind)),da);
    NLy = g2k( 2*S(:,:,v1_ind).*S(:,:,v0_ind) + sqr2*S(:,:,v2_ind).*S(:,:,v1_ind),da);
    NLz = g2k( S(:,:,v0_ind).*w1 + sqr2i*(S(:,:,v1_ind).*w2 - S(:,:,v2_ind).*w1) ,da);

    Lk = -f * Sk(:,:,u1_ind) + 1/m*iky_.*Sk(:,:,b1_ind);

    Rk(:,:,v1_ind,1) = Lk - (ikx_.*NLx + iky_.*NLy + m*NLz);
    

    % ----- RHS u2 ----- 
    NLx = g2k(2*S(:,:,u2_ind).*S(:,:,u0_ind) + sqr2i*S(:,:,u1_ind).*S(:,:,u1_ind) ,da);
    NLy = g2k(  S(:,:,v0_ind).*S(:,:,u2_ind) + S(:,:,v2_ind).*S(:,:,u0_ind) + ...
           sqr2i*S(:,:,v1_ind).*S(:,:,u1_ind) ,da);
    NLz = g2k( 2*S(:,:,u0_ind).*w2 + sqr2*S(:,:,u1_ind).*w1 ,da);

    Lk = f * Sk(:,:,v2_ind) +  (1/(2*m))*ikx_.*Sk(:,:,b2_ind);

    Rk(:,:,u2_ind,1) = Lk - (ikx_.*NLx + iky_.*NLy + m*NLz);
    
    % ----- RHS v2 ----- 
    NLx = g2k( S(:,:,u0_ind).*S(:,:,v2_ind) + S(:,:,u2_ind).*S(:,:,v0_ind)...
               + sqr2i*S(:,:,u1_ind).*S(:,:,v1_ind) ,da);
    NLy = g2k( 2*S(:,:,v2_ind).*S(:,:,v0_ind)...
               + sqr2i*S(:,:,v1_ind).*S(:,:,v1_ind) ,da);
    NLz = g2k( 2*S(:,:,v0_ind).*w2 + sqr2*S(:,:,v1_ind).*w1 ,da);

    Lk = - f * Sk(:,:,u2_ind) + (1/(2*m))*iky_.*Sk(:,:,b2_ind);
    
    Rk(:,:,v2_ind,1) = Lk - 1*(ikx_.*NLx + iky_.*NLy + m*NLz);
    
    
    % ----- RHS b1 -----
    NLx = g2k( S(:,:,u0_ind).*S(:,:,b1_ind) + ...
          sqr2i*S(:,:,u1_ind).*S(:,:,b2_ind) - sqr2i*S(:,:,u2_ind).*S(:,:,b1_ind) ,da);
    NLy = g2k( S(:,:,v0_ind).*S(:,:,b1_ind) + ...
          sqr2i*S(:,:,v1_ind).*S(:,:,b2_ind) - sqr2i*S(:,:,v2_ind).*S(:,:,b1_ind) ,da);
    NLz = -sqr2i*g2k( S(:,:,b2_ind).*w1 + S(:,:,b1_ind).*w2 ,da);

    Lk = - w1k;

    Rk(:,:,b1_ind,1) = N^2 * Lk - 1*(ikx_.*NLx + iky_.*NLy + m*NLz);
    
    
    % ----- RHS b2 ----- 
    NLx = g2k( S(:,:,u0_ind).*S(:,:,b2_ind) +  sqr2i*S(:,:,u1_ind).*S(:,:,b1_ind) ,da);
    NLy = g2k( S(:,:,v0_ind).*S(:,:,b2_ind) +  sqr2i*S(:,:,v1_ind).*S(:,:,b1_ind) ,da);
    NLz =  sqr2*g2k( S(:,:,b1_ind).*w1 ,da);
    
    Lk = -w2k;
    
    Rk(:,:,b2_ind,1) = N^2 * Lk -  (ikx_.*NLx + iky_.*NLy + m*NLz);
   
    % ----- Max velocity for CFL check ----- 
    vel_ind = [u0_ind v0_ind u1_ind u2_ind v1_ind v2_ind];
    %Umax = max([max(max(max(abs(S(:,:,vel_ind))))) max(max(abs(u0))) max(max(abs(v0))) Cmax]);
    Umax = max([max(max(max(abs(S(:,:,vel_ind))))) Cmax]);
    
    % Testing    
    Rk(:,:,kill_list,:) = 0;
    Diag2.time(counter) = t;
    Diag2.ke0(counter) = sum(sum(abs(Sk(:,:,u0_ind)).^2 + abs(Sk(:,:,v0_ind)).^2 )) ...
                        - 0.5*abs(Sk(1,1,u0_ind))^2 - 0.5 * abs(Sk(1,1,v0_ind))^2;
    Diag2.ke1(counter) = sum(sum(abs(Sk(:,:,u1_ind)).^2 + abs(Sk(:,:,v1_ind)).^2 )) ...
                                  -0.5*abs(Sk(1,1,u1_ind))^2 - 0.5*abs(Sk(1,1,v1_ind))^2;
    Diag2.ke2(counter) = sum(sum(abs(Sk(:,:,u2_ind)).^2 + abs(Sk(:,:,v2_ind)).^2))...
                                  -0.5*abs(Sk(1,1,u2_ind))^2 - 0.5*abs(Sk(1,1,v2_ind))^2;
    Diag2.kev2(counter) = sum(sum(abs(Sk(:,:,v2_ind)).^2)) - 0.5*abs(Sk(1,1,v2_ind))^2;
    Diag2.keu2(counter) = sum(sum(abs(Sk(:,:,u2_ind)).^2)) - 0.5*abs(Sk(1,1,u2_ind))^2;    
    Diag2.pe1(counter) = 1/N^2*sum(sum(abs(Sk(:,:,b1_ind)).^2))  -1/N^2 * 0.5*abs(Sk(1,1,b1_ind))^2; 
    Diag2.pe2(counter) = 1/N^2*sum(sum(abs(Sk(:,:,b2_ind)).^2)) - 1/N^2 * 0.5*abs(Sk(1,1,b2_ind))^2;
    Diag2.ke(counter)  = Diag2.ke0(counter) + Diag2.ke1(counter) + Diag2.ke2(counter);
    Diag2.pe(counter) = Diag2.pe1(counter) + Diag2.pe2(counter);
    
    Diag2.dke0dt(counter) = -mean(mean( S(:,:,u0_ind).* ( S(:,:,u1_ind).*k2g(ikx_.* Sk(:,:,u1_ind),da) +...
                                              S(:,:,v1_ind).*k2g(iky_.* Sk(:,:,u1_ind),da) ) +...
                            S(:,:,v0_ind).* ( S(:,:,u1_ind).*k2g(ikx_.* Sk(:,:,v1_ind),da) + ...
                                              S(:,:,v1_ind).*k2g(iky_.* Sk(:,:,v1_ind),da) ) ));
    Diag2.dke0dt(counter) =  Diag2.dke0dt(counter) ...
                            + m * mean(mean(w1.*(S(:,:,u0_ind).*S(:,:,u1_ind) + S(:,:,v0_ind).*S(:,:,v1_ind))));
    Diag2.bflux1(counter) = -mean(mean(w1.*S(:,:,b1_ind)));
    
    % Potential vorticity. It's computed here to avoid computing w again
   % if (mod(counter-1,frame)==0)
        
    %    v0x = k2g(ikx_ .* Sk(:,:,v0_ind), noda);
    %    v1x = k2g(ikx_ .* Sk(:,:,v1_ind), noda);
    %    b1x = k2g(ikx_ .* Sk(:,:,b1_ind), noda);
    %    v2x = k2g(ikx_ .* Sk(:,:,v2_ind), noda);
    %    b2x = k2g(ikx_ .* Sk(:,:,b2_ind), noda);
    %    w1x = k2g(ikx_ .* w1k, noda);
    %    w2x = k2g(ikx_ .* w2k, noda);
        
    %    u0y = k2g(iky_ .* Sk(:,:,u0_ind), noda);
    %    u1y = k2g(iky_ .* Sk(:,:,u1_ind), noda);
    %    b1y = k2g(iky_ .* Sk(:,:,b1_ind), noda);
    %    u2y = k2g(iky_ .* Sk(:,:,u2_ind), noda);
    %    b2y = k2g(iky_ .* Sk(:,:,b2_ind), noda);
    %    w1y = k2g(iky_ .* w1k, noda);
    %    w2y = k2g(iky_ .* w2k, noda);
        
    %    u1 = k2g(Sk(:,:,u1_ind), noda);
    %    v1 = k2g(Sk(:,:,v1_ind), noda);
    %    b1 = k2g(Sk(:,:,b1_ind), noda);
    %    u2 = k2g(Sk(:,:,u2_ind), noda);
    %    v2 = k2g(Sk(:,:,v2_ind), noda);
    %    b2 = k2g(Sk(:,:,b2_ind), noda);
        
    %    q(:,:,1) = N^2 * (v0x - u0y) +  ( b1x .* (w1y +   m * v1) - b1y .* (  m * u1 + w1x) +   m * b1 .* (v1x - u1y) + ...
     %                                     b2x .* (w2y + 2*m * v2) - b2y .* (2*m * u2 + w2x) + 2*m * b2 .* (v2x - u2y));
    %    q(:,:,2) = f * m * b1 + N^2 * (v1x - u1y) + m * b1 .* (v0x - u0y) ...
     %                  +  1/2 * sqrt(2) * (         b2x .* (w1y + m * v1) +    b1x .* (w2y + 2 * m * v2)...
     %                  -       -b2y .* (w1x + m * u1) -    b1y .* (w2x + 2 * m * u2)...
     %                  + 2 * m * b2  .* (v1x - u1y)    + m * b1 .* (v2x - u2y) );
     %   q(:,:,3) = f * 2 * m * b2 + N^2 * (v2x - u2y) + 2 *m *b2 .* (v0x - u0y)...
     %                  + 1/2 * sqrt(2) * ( -b1x .* (w1y + m * v1) + b1y .* (w1x + m * u1) + m * b1 .* (v1x - u1y) );
   % end
    
          
end
    
%-------------------------------------------------------------------

function [fg] = k2g(fk,dealias)
    
normfac = ( size(fk,1) )^2;
if dealias, normfac = (3/2)^2*normfac; end
fg = normfac*real(ifft(ifft(sympad(fk,dealias),[],1),[],2));

end 

%-------------------------------------------------------------------
    
function [fk] = g2k(fg,dealias)
    
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

fk(n+1:end,1,:) = 0;  % this part given by C-S, so is redundant
% Modified by WD on  NOV 5, 2019
% remove nx/2+1 modes
fk(n+1,:) = 0;
end

%-------------------------------------------------------------------

function [hkb] = sympad(hk,dealias)

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

%-------------------------------------------------------------------  

end % end of entire function

