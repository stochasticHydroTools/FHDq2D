function [U, V] = RandomVelocity2(NX, NY, LX, LY, KX, KY, dx, dy, epsilon, K, C1, C2)
   % Generates a quasi-2D random velocity field
   % epsilon = k_B*T/eta 
   % Generate velocity spectrum in Fourier space:
   % K = sqrt(Ksquare); 
   % C1, C2 kernel functions
   
   % Generate two scalar white noise fields: in Fourier space:
   Zk1 = Zgen(NY,NX);
   Zk2 = Zgen(NY,NX);

   % Now compute velocity in Fourier space:
   Uhat = sqrt(2*epsilon)*( sqrt(C1).*(KX).*Zk1 + sqrt(C2).*(-KY).*Zk2 )./(K.^(3/2));
   Vhat = sqrt(2*epsilon)*( sqrt(C1).*(KY).*Zk1 + sqrt(C2).*(KX).*Zk2 )./(K.^(3/2));
   
   % Enforce real-valued velocity by imposing conjugate symmetry:
   Uhat = ConjSym(Uhat,0);
   Vhat = ConjSym(Vhat,0);
   
   % fourier transform and scale to convert to real space
   U = fs2real2d(Uhat,dx,dy);
   V = fs2real2d(Vhat,dx,dy);  

   %%% uncomment this and you can have the velocity field plot
   if(0)
     figure(10);clf
     quiver(U,V)
     title('Random velocity')
   end

end  %% for function def
