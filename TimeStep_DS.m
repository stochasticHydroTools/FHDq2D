function [output] = TimeStep_DS(chat,dt,dx,dy,NX,NY,LX,LY,KX,KY,iKX,iKY,Ksquare,...
    c_whole_hat,c,c_background_modified,velMode,chi,NoVel,NoModified,epsilon,avgr,alpha,NoConvl,U,V,RXhat,RYhat,filter_type,filter)
% Takes a step for the stochastic advection-diffusion 
% Uses Crank-Nicolson for the diffusion but forward Euler for the random advection, noise and convolution
% epsilon = k_B*T/eta
% avgr is value for average gradient(const) of c, used for linearized FHD computations

source = zeros(NY,NX); % Right-hand side for implicit diffusion solve 

% Advection by random velocity div(w*c)
if(NoModified && ~NoVel) %% add in velocity term for full eqn
    source = source - (iKX.*real2fs2d(U.*c, dx, dy) + ...
                       iKY.*real2fs2d(V.*c, dx, dy)); 
end   

% For linearized FHD we linearize advection terms:
if(~NoModified) %% add in velocity term for modified eqn 
    source = source - real2fs2d(V*avgr, dx, dy);
    source = source - (real2fs2d(alpha * c_background_modified.*U, dx, dy).*iKX+...
                       real2fs2d(alpha * c_background_modified.*V, dx, dy).*iKY);
end



% Add in stochastic mass flux sqrt(2*chi*c)*white_noise
if(alpha)
   % Generate white noise in real space:
   WX = randn(NY,NX)/(sqrt(dx*dy)); WY = randn(NY,NX)/(sqrt(dx*dy));
   
   % Scale by sqrt(2*chi*c):
   if(NoModified) % Nonlinear FHD
     amplitude=sqrt(2*alpha*max(0,c)*chi);
     if(filter_type>0) 
         WX = fs2real2d(real2fs2d(WX, dx, dy).*filter, dx, dy);
         WY = fs2real2d(real2fs2d(WY, dx, dy).*filter, dx, dy);
     end
     
   elseif(~NoModified) % Linearized FHD
     amplitude=sqrt(2*alpha*c_background_modified*chi);
   end
   WX = WX.*amplitude; WY = WY.*amplitude;

   % Take it to Fourier space to compute div(stochastic_flux):
   WXhat = real2fs2d(WX, dx, dy); WYhat = real2fs2d(WY, dx, dy);
   source = source + iKX.*WXhat + iKY.*WYhat;
end

% Terms added so far are all white in time:
source = source/sqrt(dt); % scaling for random noise and random velocity


% For Quasi2D, we have an additional nonlinear term to compute
if(velMode==1) % add in convolution term
    RCXhat = RXhat.*c_whole_hat; RCYhat = RYhat.*c_whole_hat;
    VRCX = fs2real2d(RCXhat, dx, dy); VRCY = fs2real2d(RCYhat, dx, dy); 
    if(~NoModified) % Linearized FHD:
      source = source + epsilon*(iKX.*real2fs2d(c_background_modified.*VRCX, dx, dy) + ...
                                 iKY.*real2fs2d(c_background_modified.*VRCY, dx, dy)); 
    elseif(NoModified && ~NoConvl) % Nonlinear FHD
       source = source + epsilon*(iKX.*real2fs2d(VRCX.*c, dx, dy) + ...
                                  iKY.*real2fs2d(VRCY.*c, dx, dy)); 
    end    
end

if(filter_type<0) 
    source = source.*filter;
end

% Do implicit Crank-Nicolson solve for diffusion spectrally:
output = ((1 - 0.5*dt*chi*Ksquare).*chat + source*dt)./(1 + 0.5*dt*chi*Ksquare);
