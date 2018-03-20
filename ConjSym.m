function output = ConjSym( X,x11 )
   %   the function conjugate-symmetrize a matrix for further ifft2
   %   intput X should be a square matrix, wavenumber order is
   %   0,1,...,N/2-1,-N/2,-N/2+1,..,-1 for both row and column
   %   X(1,1)=x11 should be a real number

   [NY,NX] = size(X);

   output = X;

   output(2:NY,NX/2+2:NX) = conj(output(NY:-1:2,NX/2:-1:2)); 
   output(1,NX/2+2:NX) = conj(output(1,NX/2:-1:2)); 
   output(NY/2+2:NY,1) = conj(output(NY/2:-1:2,1)); 
   output(NY/2+2:NY,NX/2+1) = conj(output(NY/2:-1:2,NX/2+1));
   output(1,1)=x11;
   
end
