function [outputX,outputY] = R_dot_grad( iKX, iKY, k, C1)
   %   This produce RnablaX, RnablaY for convolution 
   

   Rtemp = C1;
   
   outputX = Rtemp .* iKX./k;
   outputY = Rtemp .* iKY./k;


end
