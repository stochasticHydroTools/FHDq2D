function output = Zgen( NY,NX )
   % Generates a scalar white noise field in Fourier space
   % This is not as trivial as doing it in real space because it's a complex number

   output = sqrt(0.5)*randn( NY,NX )+1i*sqrt(0.5)*randn( NY,NX );

   % Ensure the result is real by treating the special mode without a conjugate partner separately:
   output(NY/2+1,1)=randn;
   output(1,NX/2+1)=randn;
   output(NY/2+1,NX/2+1)=randn;

end

