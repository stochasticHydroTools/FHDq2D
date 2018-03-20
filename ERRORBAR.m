

%%% this file is used to calculate mean and errorbar data for spectrum  
clear all;
origin_filename = 'eps1vel1conv0noise0color_Num2VM2';
load([origin_filename,'/parameters.mat'])

filename = [origin_filename,'/',num2str(NX),'eps',num2str(epsilon),'vm',num2str(velMode),'dt',num2str(dt)];

KX=(2*pi/LX)*ones(1,NY)'*(mod((1:NX)-ceil(NX/2+1),NX)-floor(NX/2));
KX_1d = KX(1,2:NX/2); %% wavenumber for x coordinate


specGR_1 = zeros(rerunNum,(NX-2)/2);specGR_2 = zeros(rerunNum,(NX-2)/2);
specGG_1 = zeros(rerunNum,(NX-2)/2);specGG_2 = zeros(rerunNum,(NX-2)/2);
specRR_1 = zeros(rerunNum,(NX-2)/2);specRR_2 = zeros(rerunNum,(NX-2)/2);

%%%% read in data and take means and std
for j = 1:rerunNum  %%% GR
    load([filename,'/specGR_record_half_1_',num2str(j),'.mat']);
    load([filename,'/specGR_record_half_2_',num2str(j),'.mat']);
    specGR_1(j,:) = -specGR_record_half_1(1:(NX-2)/2)/epsilon;
    specGR_2(j,:) = -specGR_record_half_2(1:(NX-2)/2)/epsilon;
end

for j = 1:rerunNum    %%%% GG
    load([filename,'/specGG_record_half_1_',num2str(j),'.mat']);
    load([filename,'/specGG_record_half_2_',num2str(j),'.mat']);
    specGG_1(j,:) = specGG_record_half_1(1:(NX-2)/2)/epsilon;
    specGG_2(j,:) = specGG_record_half_2(1:(NX-2)/2)/epsilon;
end

for j = 1:rerunNum    %%%% RR
    load([filename,'/specRR_record_half_1_',num2str(j),'.mat']);
    load([filename,'/specRR_record_half_2_',num2str(j),'.mat']);
    specRR_1(j,:) = specRR_record_half_1(1:(NX-2)/2)/epsilon;
    specRR_2(j,:) = specRR_record_half_2(1:(NX-2)/2)/epsilon;
end


specGR_1_mean = mean(specGR_1);specGR_2_mean = mean(specGR_2);
specGR_1_std = std(specGR_1);specGR_2_std = std(specGR_2);
specGR_1_error = 2*specGR_1_std/sqrt(rerunNum);specGR_2_error = 2*specGR_2_std/sqrt(rerunNum);

specGG_1_mean = mean(specGG_1);specGG_2_mean = mean(specGG_2);
specGG_1_std = std(specGG_1);specGG_2_std = std(specGG_2);
specGG_1_error = 2*specGG_1_std/sqrt(rerunNum);specGG_2_error = 2*specGG_2_std/sqrt(rerunNum);

specRR_1_mean = mean(specRR_1);specRR_2_mean = mean(specRR_2);
specRR_1_std = std(specRR_1);specRR_2_std = std(specRR_2);
specRR_1_error = 2*specRR_1_std/sqrt(rerunNum);specRR_2_error = 2*specRR_2_std/sqrt(rerunNum);

%%%%% save data
%%%% GR
specGR_1_data_eps1_64_anti = [KX_1d' specGR_1_mean' specGR_1_error'];
specGR_2_data_eps1_64_anti = [KX_1d' specGR_2_mean' specGR_2_error'];
save([origin_filename,'/specGR_1_data_eps',num2str(epsilon),'_',num2str(NX),'.txt'],'specGR_1_data_eps1_64_anti','-ascii');
save([origin_filename,'/specGR_2_data_eps',num2str(epsilon),'_',num2str(NX),'.txt'],'specGR_2_data_eps1_64_anti','-ascii');

%%% GG
specGG_1_data_eps1_64_anti = [KX_1d' specGG_1_mean' specGG_1_error'];
specGG_2_data_eps1_64_anti = [KX_1d' specGG_2_mean' specGG_2_error'];
save([origin_filename,'/specGG_1_data_eps',num2str(epsilon),'_',num2str(NX),'.txt'],'specGG_1_data_eps1_64_anti','-ascii');
save([origin_filename,'/specGG_2_data_eps',num2str(epsilon),'_',num2str(NX),'.txt'],'specGG_2_data_eps1_64_anti','-ascii');

% %%%% RR
specRR_1_data_eps1_64_anti = [KX_1d' specRR_1_mean' specRR_1_error'];
specRR_2_data_eps1_64_anti = [KX_1d' specRR_2_mean' specRR_2_error'];
save([origin_filename,'/specRR_1_data_eps',num2str(epsilon),'_',num2str(NX),'.txt'],'specRR_1_data_eps1_64_anti','-ascii');
save([origin_filename,'/specRR_2_data_eps',num2str(epsilon),'_',num2str(NX),'.txt'],'specRR_2_data_eps1_64_anti','-ascii');

%%%% plot first and second half data 
if(1)
   if(filter) %% anti-aliasing
    NX = 2*NX/3;
    NY = 2*NY/3;
   end 
   
    figure(1);clf
    %%% GR 
    errorbar(specGR_1_data_eps1_64_anti(1:(NX-2)/2,1),specGR_1_data_eps1_64_anti(1:(NX-2)/2,2),specGR_1_data_eps1_64_anti(1:(NX-2)/2,3));hold on
    errorbar(specGR_2_data_eps1_64_anti(1:(NX-2)/2,1),specGR_2_data_eps1_64_anti(1:(NX-2)/2,2),specGR_2_data_eps1_64_anti(1:(NX-2)/2,3));
    set(gca,'xscale','log'); set(gca,'yscale','log');    
    
    %%% GG
    errorbar(specGG_1_data_eps1_64_anti(1:(NX-2)/2,1),specGG_1_data_eps1_64_anti(1:(NX-2)/2,2),specGG_1_data_eps1_64_anti(1:(NX-2)/2,3));
    errorbar(specGG_2_data_eps1_64_anti(1:(NX-2)/2,1),specGG_2_data_eps1_64_anti(1:(NX-2)/2,2),specGG_2_data_eps1_64_anti(1:(NX-2)/2,3));
   
    %%% RR
    errorbar(specRR_1_data_eps1_64_anti(1:(NX-2)/2,1),specRR_1_data_eps1_64_anti(1:(NX-2)/2,2),specRR_1_data_eps1_64_anti(1:(NX-2)/2,3));
    errorbar(specRR_2_data_eps1_64_anti(1:(NX-2)/2,1),specRR_2_data_eps1_64_anti(1:(NX-2)/2,2),specRR_2_data_eps1_64_anti(1:(NX-2)/2,3));
     
    
    legend('GR-1','GR-2','GG-1','GG-2','RR-1','RR-2')
    if(velMode==1)%% Q2D
        load('Q2D_first.dat');
        load('Q2D_second.dat');
        errorbar(Q2D_first(1:NX/2-1,1),Q2D_first(1:NX/2-1,2),Q2D_first(1:NX/2-1,3));
        errorbar(Q2D_second(1:NX/2-1,1),Q2D_second(1:NX/2-1,2),Q2D_second(1:NX/2-1,3));
        legend('GR-1','GR-2','GG-1','GG-2','RR-1','RR-2','Q2D-1','Q2D-2')
    elseif(velMode==2)%%% T2D
        load('T2D_first.dat');
        load('T2D_second.dat');
        errorbar(T2D_first(1:NX/2-1,1),T2D_first(1:NX/2-1,2),T2D_first(1:NX/2-1,3));
        errorbar(T2D_second(1:NX/2-1,1),T2D_second(1:NX/2-1,2),T2D_second(1:NX/2-1,3));
        legend('GR-1','GR-2','GG-1','GG-2','RR-1','RR-2','T2D-1','T2D-2')
    end
end    

