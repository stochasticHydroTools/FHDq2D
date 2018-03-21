function [c1,c2] = VelKernel(hydro_kernel, k, sigma, Ksquare, kc)
% Compute the functions c_1(k) and c_2(k) based on the type of hydrodynamics
%   hydro_kernel: 1 for quasi2D, 2 for true2D, 3 for Saffman
%                -1 for quasi2D with incompressible component only(c1=0)
%   Ksquare = KX.^2 + KY.^2
%   k = sqrt(Ksquare)

K = k*sigma;

if(hydro_kernel==1)
    c1 = 0.25*(-2*K.*exp(-Ksquare*sigma^2/pi) - (pi+2*K.^2).*( erf(K/sqrt(pi)) -1 ) )/pi;
    c2 = -0.5*erf(K/sqrt(pi))+0.5;
elseif(hydro_kernel==-1)
    c1 = zeros(size(K)); 
    c2 = -0.5*erf(K/sqrt(pi))+0.5;
elseif(hydro_kernel==2)
    c1 = zeros(size(K));
    c2 = exp(-(Ksquare*sigma^2)/pi)./(k); % True2D is the same as Saffman with kc=0
elseif(hydro_kernel==3)
    c1 = zeros(size(K));
    c2 = exp(-(Ksquare*sigma^2)/pi)./(k+kc);
end
        
end

