function output_c = initial_2color(c_background, proportion, equilibrium,NX,NY)
% Generates initial conditions corresponding to a stripe (nonequilibrium) or uniform field (equilibrium)

c = zeros(NY,NX,2);

if(equilibrium)
    c(:,:,1) = c(:,:,1) + c_background/2;
    c(:,:,2) = c(:,:,2) + c_background/2;
else
    if(mod(NY-1,3)==0)
        c(:,:,1) = c(:,:,1) + c_background;
        c((NY-1)/3 + 1,:,1) = (1+proportion)/2*c_background;
        c((NY-1)*2/3 + 1,:,1) = (1+proportion)/2*c_background;
        c((NY-1)/3 + 2:(NY-1)*2/3,:,1) = proportion*c_background;
        
        c(:,:,2) = c(:,:,2) + c_background*proportion;
        c((NY-1)/3 + 1,:,2) = (1+proportion)/2*c_background;
        c((NY-1)*2/3 + 1,:,2) = (1+proportion)/2*c_background;
        c((NY-1)/3 + 2:(NY-1)*2/3,:,2) = c_background;
        
    elseif(mod(NY-1,3)==1)
        c(:,:,1) = c(:,:,1) + c_background;
        c(floor((NY-1)/3) + 1,:,1) = (5/6+proportion/6)*c_background;
        c(NY-floor((NY-1)/3),:,1) = (5/6+proportion/6)*c_background;
        c(floor((NY-1)/3) + 2:NY-floor((NY-1)/3)-1,:,1) = proportion*c_background;
        
        c(:,:,2) = c(:,:,2) + c_background*proportion;
        c(floor((NY-1)/3) + 1,:,2) = (1/6+proportion*5/6)*c_background;
        c(NY-floor((NY-1)/3),:,2) = (1/6+proportion*5/6)*c_background;
        c(floor((NY-1)/3) + 2:NY-floor((NY-1)/3)-1,:,2) = c_background;
    elseif(mod(NY-1,3)==2)
        c(:,:,1) = c(:,:,1) + c_background;
        c(floor((NY-1)/3) + 2,:,1) = (1/6+proportion*5/6)*c_background;
        c(NY-floor((NY-1)/3)-1,:,1) = (1/6+proportion*5/6)*c_background;
        c(floor((NY-1)/3) + 3:NY-floor((NY-1)/3)-2,:,1) = proportion*c_background;
        
        c(:,:,2) = c(:,:,2) + c_background*proportion;
        c(floor((NY-1)/3) + 2,:,2) = (5/6+proportion/6)*c_background;
        c(NY-floor((NY-1)/3)-1,:,2) = (5/6+proportion/6)*c_background;
        c(floor((NY-1)/3) + 3:NY-floor((NY-1)/3)-2,:,2) = c_background;
    end      
end

output_c = c;

end

