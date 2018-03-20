function  output = fs2real2d( fhat, dx, dy )
% Converts a field from Fourier to real space with the correct scaling

 NY = size(fhat,1);
 NX = size(fhat,2);
 output = ifft2(fhat)*sqrt(NY*NX)/sqrt(dx*dy);

end

