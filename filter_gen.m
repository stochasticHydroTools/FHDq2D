function [ output ] = filter_gen( NY, NX, filter_type )
%  filter_type: 0 for no filter, 1 for 1/2 filter, 2 for 2/3 filter
%  Assuming no fftshift used
   output = ones(NY,NX);
   
   if(filter_type == 1)
       output(1+NY/2-1/4*NY:1+NY/2+1/4*NY,:) = 0;
       output(:,1+NX/2-1/4*NX:1+NX/2+1/4*NX) = 0;
   elseif(filter_type == 2)
       output(1+NY/2-1/6*NY:1+NY/2+1/6*NY,:) = 0;
       output(:,1+NX/2-1/6*NX:1+NX/2+1/6*NX) = 0;
   end
      
end

