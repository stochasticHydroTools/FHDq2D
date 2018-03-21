clear all; 

%%%%% parameters %%%%%%%%%%%
%---------------------------
NoVel = 0; % Exclude random advection in full eqn
NoConvl =1; %% Exclude convolution term in full eqn
alpha = 0; % parameter for modified equation, also impose noise term in full equation if not 0
NoModified = 1; %% indicate we are using linearized eqn or not
avgr = 0; % avgr is average gradient of c for linearized/modified equations
%%% general parameters %%%%%%%%%
velMode = -1; %% 1 for quasi2D and 2 for true2D velocity type, 3 for Saffman 
             %%% -1 for quasi2D with incompressible component only 
%kc = 0; %%% const for Saffman 
%kc = 0.0224193; kc = 0.0448385; kc = 0.0896771; kc = 0.49995;
kc = 0.179354;

epsilon = 1e0; % epsilon = k_B*T/eta

%chi = 0.3993036270; % Diffusion coefficient for True2D
chi = 0.857/(6*pi); % Diffusion coefficient for Quasi2D
%chi = 0.399315; % Diffusion coefficient for Saffman kc=0
%chi = 0.292203; % Diffusion coefficient for Saffman kc=0.022
%chi = 0.248638; % Diffusion coefficient for Saffman kc=0.044
%chi = 0.200213; % Diffusion coefficient for Saffman kc=0.089
%chi = 0.148534; % Diffusion coefficient for Saffman kc=0.179
%chi = 0.0686706; % Diffusion coefficient for Saffman kc=0.499

c_background = 0.3183088885; 
uniform = 0;  %% Put a bump or just do a homogeneous initial 
proportion = 0; %% proportion of bump if exists

%NX= 32; NY = 32; %% number of grids
%NX = 64; NY = 64;
%NX = 48; NY = 48;
NX = 96; NY = 96;
%NX = 128; NY = 128;
%NX = 192; NY = 192;
%NX = 256; NY = 256;

LX = 560.5; LY = 560.5; %% length of area
dx = LX/(NX-1); dy = LY/(NY-1);
%dx = LX/(NX); dy = LY/(NY);
%dt = 10; %% CFL~0.5
CFL = 0.5; dt =  CFL*dx^2/chi;  %% setting dt with CFL number

T = 1.505e5; %Quasi2D
%T = 20000; % True2D
%T = 17530; % kc = 0 Saffman
%T = 23955.9; % kc = 0.022 Saffman
%T = 28153.4; % kc = 0.044 Saffman
%T = 34962.8; % kc = 0.089 Saffman
%T = 47127.3; % kc = 0.179 Saffman
%T = 101936; % kc = 0.499 Saffman
MaxIter = floor(T/dt);  %% number of iterations
rerunNum = 16 ;  %% times to rerun the simulation(for stage test)

color_Num = 2; % number of colors
filter = 1; % use 2/3 filter or antialising?
clims = [0 c_background];
sigma = 1; % Hydrodynamic radius of particles
if(0)  %% set sigma by approximation(8) in notes from chi
    if(velMode==1) %% quasi2D
       sigma = epsilon/(6*pi*chi);
    else  %% true2D
       sigma = LX/(3.71*exp(4*pi*chi/epsilon));
    end
end
filename = ''; % prefix for file names

% Control the run and what is output:
% when gap is 0, collection of that data is skipped
real2d_gap = floor(MaxIter/50); %% real 2d evolution movie
real_spec_2dMean_gap =  1; %% snapshot of real 2d mean and 2d spectrum
spectrum1d_gap = 1; %% snapshot of 1d spectrum along x
CYmean_gap = 0;  %% snapshot of c(y), mean of c on x  

% ------------------------------------
%%% initialize statistics variables

if(real_spec_2dMean_gap)
    cMean_record = zeros(MaxIter/real_spec_2dMean_gap+1,color_Num+1); %% save mean of real 2d
    spec2dmean = zeros(MaxIter/real_spec_2dMean_gap+1,color_Num+1); %% save mean of 2d spectrum
end

if(spectrum1d_gap)
    spec1dMean_record = zeros(MaxIter/spectrum1d_gap+1); %% save 1d spectrum mean
    cVar_record =  zeros(MaxIter/spectrum1d_gap+1); %% save real variance
    specX_record = zeros((NX-2)/2,MaxIter/spectrum1d_gap+1); %% save 1d spectrum of density along x
    if(color_Num==2)
        specGG_record = zeros((NX-2)/2,MaxIter/spectrum1d_gap+1); %% save 1d green spectrum along x
        specRR_record = zeros((NX-2)/2,MaxIter/spectrum1d_gap+1); %% save 1d red spectrum along x
        specGR_record = zeros((NX-2)/2,MaxIter/spectrum1d_gap+1); %% save 1d green-red spectrum along x
        specGR_record_half_1 = zeros((NX-2)/2,1);
        specGR_record_half_2 = zeros((NX-2)/2,1);
    end
end
if(CYmean_gap)
   CYmean_record = zeros(NY,MaxIter/CYmean_gap+1); %% save c(y), mean of c on x
   if(color_Num==2)
       CYmean_record_red = zeros(NY,MaxIter/CYmean_gap+1); %% save c(y), mean of c on x
       CYmean_record_green = zeros(NY,MaxIter/CYmean_gap+1); %% save c(y), mean of c on x
   end
end

%%% create folder for data saving
if(NoModified)
    filename = [filename,'eps',num2str(epsilon),'vel',num2str(1-NoVel),'conv',num2str(1-NoConvl),'noise',num2str(alpha),'color_Num',num2str(color_Num),'VM',num2str(velMode)];
else
    filename = [filename,'Modified_eps',num2str(epsilon),'vel',num2str(1-NoVel),'conv',num2str(1-NoConvl),'noise',num2str(alpha),'color_Num',num2str(color_Num),'VM',num2str(velMode)];
end
mkdir(filename);


% ------------------------------------
%%% Pre-compute quantities reused inside time loop many times

N_particles_per_cell = c_background*dx*dy;

KX=(2*pi/LX)*ones(1,NY)'*(mod((1:NX)-ceil(NX/2+1),NX)-floor(NX/2)); % matrix of wavenumbers in x direction 
KY=(2*pi/LY)*(mod((1:NY)'-ceil(NY/2+1),NY)-floor(NY/2))*ones(1,NX); % matrix of wavenumbers in y direction 
Ksquare = KX.^2+KY.^2;
K = sqrt(Ksquare); %k
K(1,1)=1; % Avoid division by zero


% For spectral differentiation we need to remove the unmatched wavenumber:
iKX=1i*KX;
iKY=1i*KY;
%%% remove the 'largest' mode to make sure ifft2 gives real output
if(mod(NX,2)==0) % Even grid size
   iKX(:,NX/2+1)=0;  
end
if(mod(NY,2)==0) % Even grid size
   iKY(NY/2+1,:)=0;
end


[C1, C2] = VelKernel(velMode, K, sigma, Ksquare, kc);

[RXhat,RYhat] = R_dot_grad( iKX,iKY,  K, C1); %% convolution in fourier


% ------------------------------------
%%% deterministic part of initial condition (reused for each independent run)

c = zeros(NY,NX,color_Num);
if(color_Num==1)
if(uniform)
    c = c + c_background;
else
    c = c + c_background*proportion;
    c (floor(NY/3):floor(2*(NY/3)),:) = c_background;
end
elseif(color_Num==2) %% c(:,:,1) is green and c(:,:,2) is red
if(~NoModified)
c(:,:,1) = c(:,:,1) + c_background/2;
c(:,:,2) = c(:,:,2) + c_background/2;
else
c = initial_2color(c_background, proportion, 0,NX,NY);
end
end

c_noNoise = c;
if(color_Num==1)
    c_background_modified = c_background;
else
    c_background_modified = c_background/2;
end

% ------------------------------------
%%% Begin repeat loop for doing multiple simulations and averaging over them

mkdir([filename,'/',num2str(NX),'eps',num2str(epsilon),'vm',num2str(velMode),'dt',num2str(dt)]);
reruntimer = 0;
Movie(floor(MaxIter/real2d_gap)+1,color_Num+1) = struct('cdata',[],'colormap',[]);
while reruntimer < rerunNum  %% rerun from the initial several times

t=0; iter = 0;  
c = c_noNoise;
chat = zeros(size(c));
specGR_record_half_1 = zeros((NX-2)/2,1); 
specGR_record_half_2 = zeros((NX-2)/2,1);
specGG_record_half_1 = zeros((NX-2)/2,1); 
specGG_record_half_2 = zeros((NX-2)/2,1);
specRR_record_half_1 = zeros((NX-2)/2,1); 
specRR_record_half_2 = zeros((NX-2)/2,1);

if(alpha~=0)  %% impose noise on initial condition
  for j = 1:color_Num
    c(:,:,j) = c(:,:,j) + randn(NY,NX).*sqrt(alpha*c(:,:,j)/(dx*dy));   
  end
end

for j = 1:color_Num
  chat(:,:,j) = real2fs2d(c(:,:,j), dx, dy);  
end

% --------------------------------------------
% Time loop iteration

while iter <= MaxIter

   for j = 1:color_Num
     c(:,:,j) = fs2real2d(chat(:,:,j), dx, dy); 
%      imagpart = max(max(imag(c)))  %% to test
%      minC = min(c(:)) %% to test

   end
   c_whole = sum(c,3);

   % ----------------------
   % HydroGrid-type spectrum analysis:
   if(real2d_gap && mod(iter,real2d_gap)==0 && reruntimer==0) %% only save movie for the first run
       Movie = Snapshot_real2d(c, c_whole, clims, iter, dt, real2d_gap, color_Num, Movie);
   end
  
   if(real_spec_2dMean_gap && mod(iter,real_spec_2dMean_gap)==0)  %% snapshot of real 2d mean and 2d spectrum
       [cMean_record, spec2dmean] = real_spec_2dmean(c, chat, iter, real_spec_2dMean_gap, NX, NY, color_Num, cMean_record, spec2dmean);
   end
   
   if(spectrum1d_gap && mod(iter,spectrum1d_gap)==0)  %% snapshot of 1d spectrum along x, and variance of c
     sumin = abs(real2fs2d(c_whole,dx,dy)).^2;
     spec1dMean_record(iter/spectrum1d_gap+1) = spec1dMean_record(iter/spectrum1d_gap+1) + mean(sumin(1,2:NX/2));
     specX_record(:,iter/spectrum1d_gap+1) = specX_record(:,iter/spectrum1d_gap+1) + sumin(1,2:NX/2)';
     cVar_record(iter/spectrum1d_gap+1) = cVar_record(iter/spectrum1d_gap+1) + var(c_whole(:)); 
     
     if(color_Num==2)
       if(iter<=MaxIter/2 && iter>=1)  %% first half of iterations
           [specGR_record, specGG_record, specRR_record,...
           specGR_record_half_1, specGG_record_half_1, specRR_record_half_1] = ...
               Spec1d_coeff(c, dx, dy, iter, NX, NY, spectrum1d_gap, specGR_record, specGG_record, specRR_record,...
               specGR_record_half_1, specGG_record_half_1, specRR_record_half_1); 
       elseif(iter>=MaxIter/2)         %% second half of iterations
           [specGR_record, specGG_record, specRR_record,...
           specGR_record_half_2, specGG_record_half_2, specRR_record_half_2] = ...
           Spec1d_coeff(c, dx, dy, iter, NX, NY, spectrum1d_gap, specGR_record, specGG_record, specRR_record,...
               specGR_record_half_2, specGG_record_half_2, specRR_record_half_2);  
       end
     end
   end
   
   if(CYmean_gap && mod(iter,CYmean_gap)==0 )  %% snapshot of c(y), mean of c on x  
     CYmean_record(:,iter/CYmean_gap+1) = CYmean_record(:,iter/CYmean_gap+1) + mean(c_whole,2); 
     if(color_Num==2)
       CYmean_record_red(:,iter/CYmean_gap+1) = CYmean_record_red(:,iter/CYmean_gap+1) + mean(c(:,:,2),2); 
       CYmean_record_green(:,iter/CYmean_gap+1) = CYmean_record_green(:,iter/CYmean_gap+1) + mean(c(:,:,1),2); 
     end
   end    
    
   
   %---------------------------
   % Update variables to time t+dt:
   
   c_whole_hat = sum(chat,3);
   [U, V] = RandomVelocity2(NX, NY, LX, LY, KX, KY, dx, dy, epsilon, K, C1, C2);  

   % This loop changes each color in turn, but only c_whole/c_hole_hat are used to update each color
   % and this does not change throughout the time step
   % Therefore this code does an Euler-Maruyuama update of chat
   for j = 1:color_Num
     chat(:,:,j) = TimeStep_DS(chat(:,:,j),dt,dx,dy,NX,NY,LX,LY,KX,KY,iKX,iKY,Ksquare,c_whole_hat,c(:,:,j),...
       c_background_modified,velMode,chi,NoVel,NoModified,epsilon,avgr*(-1)^(j+1),alpha,NoConvl,U,V,RXhat,RYhat,filter);
   end


   iter = iter+1;
   t = t+dt;

end

%----------------------------------------------
% Save statistics from run:
reruntimer = reruntimer+1;
%%% save the first half and second half of spectrum data before next run
if(color_Num >= 1)
    specGR_record_half_1 = specGR_record_half_1/(MaxIter/2); specGR_record_half_2 = specGR_record_half_2/(MaxIter/2);
    specGG_record_half_1 = specGG_record_half_1/(MaxIter/2); specGG_record_half_2 = specGG_record_half_2/(MaxIter/2);
    specRR_record_half_1 = specRR_record_half_1/(MaxIter/2); specRR_record_half_2 = specRR_record_half_2/(MaxIter/2);
    save([filename,'/',num2str(NX),'eps',num2str(epsilon),'vm',num2str(velMode),'dt',num2str(dt),'/specGR_record_half_1_',num2str(reruntimer),'.mat'],'specGR_record_half_1');
    save([filename,'/',num2str(NX),'eps',num2str(epsilon),'vm',num2str(velMode),'dt',num2str(dt),'/specGR_record_half_2_',num2str(reruntimer),'.mat'],'specGR_record_half_2');
    save([filename,'/',num2str(NX),'eps',num2str(epsilon),'vm',num2str(velMode),'dt',num2str(dt),'/specGG_record_half_1_',num2str(reruntimer),'.mat'],'specGG_record_half_1');
    save([filename,'/',num2str(NX),'eps',num2str(epsilon),'vm',num2str(velMode),'dt',num2str(dt),'/specGG_record_half_2_',num2str(reruntimer),'.mat'],'specGG_record_half_2');
    save([filename,'/',num2str(NX),'eps',num2str(epsilon),'vm',num2str(velMode),'dt',num2str(dt),'/specRR_record_half_1_',num2str(reruntimer),'.mat'],'specRR_record_half_1');
    save([filename,'/',num2str(NX),'eps',num2str(epsilon),'vm',num2str(velMode),'dt',num2str(dt),'/specRR_record_half_2_',num2str(reruntimer),'.mat'],'specRR_record_half_2');
end

end

% -----------------------------------------------
%%% Save results from averaging over multiple independent samples:

if(real2d_gap)  %%% save movies
  for j = 1:color_Num
     video = VideoWriter([filename,'/real2d_',num2str(j),'.avi']); video.FrameRate=1;
     open(video)
     writeVideo(video,Movie(:,j));
     close(video)
  end
  video = VideoWriter([filename,'/real2d_whole.avi']); video.FrameRate=1;
  open(video)
  writeVideo(video,Movie(:,color_Num+1));
  close(video)
end

if(real_spec_2dMean_gap)
    cMean_record = cMean_record/rerunNum; %% save mean of real 2d
    save([filename,'/cMean.mat'],'cMean_record');
    spec2dmean = spec2dmean/rerunNum; %% save mean of 2d spectrum
    save([filename,'/spec2dmean.mat'],'spec2dmean');
end


if(spectrum1d_gap)
    specX_record = specX_record/rerunNum; %% save 1d spectrum along x
    spec1dMean_record = spec1dMean_record/rerunNum; %% save 1d spectrum mean
    cVar_record =  cVar_record/rerunNum; %% save real variance
    if(color_Num==1)
       save([filename,'/spectrum1d.mat'],'specX_record','spec1dMean_record','cVar_record');
    else
       specGG_record = specGG_record/rerunNum;  
       specRR_record = specRR_record/rerunNum;  
       specGR_record = specGR_record/rerunNum;  
       save([filename,'/spectrum1d.mat'],'specX_record','spec1dMean_record','cVar_record','specGG_record','specGR_record','specRR_record');
    end
end
if(CYmean_gap)
   CYmean_record = CYmean_record/rerunNum; %% save c(y), mean of c on x  
   if(color_Num==1)
     save([filename,'/CYmean.mat'],'CYmean_record');
   elseif(color_Num==2)
     CYmean_record_red = CYmean_record_red/rerunNum;
     CYmean_record_green = CYmean_record_green/rerunNum;
     save([filename,'/CYmean.mat'],'CYmean_record','CYmean_record_red','CYmean_record_green');
   end
end


%%% save parameters for reference
save([filename,'/parameters.mat'],'sigma','chi','c_background','dx','dy','dt','LX','LY','NoVel','NoConvl',...
    'NoModified','velMode','epsilon','alpha','avgr','uniform','proportion','NX','NY','LX','LY','MaxIter','rerunNum',...
    'real_spec_2dMean_gap','CYmean_gap','spectrum1d_gap','color_Num','filter','c_noNoise');

%%% save parameters in a txt file for easy reading
fileID = fopen([filename,'/parameters.txt'],'w');
fprintf(fileID,'sigma = %f, chi = %12.8f\n c_background = %12.8f \n',sigma,chi,c_background);
fprintf(fileID,'color_Num = %d, uniform = %d, proportion = %f',color_Num,uniform,proportion);
fprintf(fileID,'NX=NY = %d, LX=LY = %f, dx=dy = %f, dt = %f \n',NX,LX,dx,dt);
fprintf(fileID,'MaxIter = %d, rerunNum = %d \n',MaxIter,rerunNum);
fprintf(fileID,'NoVel = %d, NoConvl = %d, NoModified = %d, filter = %d \n',NoVel,NoConvl,NoModified,filter);
if(~NoVel)
fprintf(fileID,'velMode = %d, epsilon =  %f, avgr = %f, alpha = %f \n',velMode,epsilon,avgr,alpha);
end
fclose(fileID);
