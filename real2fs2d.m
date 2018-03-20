function output= real2fs2d( f, dx, dy )
% Converts a field from real to Fourier space with the correct scaling

 NY = size(f,1);
 NX = size(f,2);
 output = fft2(f)/sqrt(NY*NX)*sqrt(dx*dy);

end

