
clear all; clc
filename = 'eps1vel1conv0noise0color_Num2VM2';%'eps0.1vel1conv1noise0color_Num2VM1_T20000'; %% path to the directory
load([filename,'/parameters.mat'])
CFL = chi*dt/dx^2;
N_particles_per_cell = c_background*dx*dy;


%%%%%% 1 to turn on the plot and 0 to turn off that plot 
real_spec_2dMean_plot = 1; %% real and spectrum 2d mean
spectrum1d_plot = 1; %% cVar, 1d spec mean plot and spectrum coefficient plot
CYmean_plot = 0; %% 
compare = 1; %% compare two CYmean from two results
filename2 = 'eps0.01vel0conv0noise0color_Num2VM1'; %% path to the result used to compare if needed




if(real_spec_2dMean_gap && real_spec_2dMean_plot)   %%% real 2d mean plot
load([filename,'/cMean.mat'])
load([filename,'/spec2dmean.mat'])
for j = 1:color_Num
 figure(j);clf
 plot(real_spec_2dMean_gap*((1:length(cMean_record(:,j)))-1),cMean_record(:,j)); hold on;
 plot(real_spec_2dMean_gap*((1:length(cMean_record(:,j)))-1),ones(1,length(cMean_record(:,j)))*mean(mean(c_noNoise(:,:,j))),'r-'); 
 plot(real_spec_2dMean_gap*((1:length(spec2dmean(:,j)))-1),spec2dmean(:,j));
 hold off; legend('real 2d mean','expected','spectrum 2d mean');xlabel('iter');
 if(color_Num==2 && j==1)
     title('mean of concentration green')
     saveas(gcf,[filename,'/real_spec_2dMean_G.png'])
 elseif(color_Num==2 && j==2)
     title('mean of concentration red')
     saveas(gcf,[filename,'/real_spec_2dMean_R.png'])
 else
     title('mean of density concentration')
     saveas(gcf,[filename,'/real_spec_2dMean_whole.png'])

 end
end
end



if(spectrum1d_gap  && spectrum1d_plot)   
load([filename,'/spectrum1d.mat'])
spectrum1d_color = 1;
figure(5);clf %% plot cVar and 1d spec mean
plot(spectrum1d_gap*((1:length(cVar_record(:,spectrum1d_color)))-1),cVar_record(:,spectrum1d_color)*(dx*dy));
% hold on; 
% plot(spectrum1d_gap*((1:length(spec1dMean_record(:,spectrum1d_color)))-1),spec1dMean_record(:,spectrum1d_color));
% hold off; legend('real 2d variance','1d spectrum mean');xlabel('iter');
% saveas(gcf,[filename,'/spectrum1d.png'])

figure(6);clf %% plot coefficient of 1d spec along x
KX=(2*pi/LX)*ones(1,NY)'*(mod((1:NX)-ceil(NX/2+1),NX)-floor(NX/2)); 
KX_1d = KX(1,2:NX/2); %% wavenumber for x coordinate
legendseries = [];
loop_gap = MaxIter/2;
for loopindex = 1:MaxIter/loop_gap
    figure(6);
    loglog(KX_1d,sum(specX_record(:,(loopindex-1)*loop_gap+1:loopindex*loop_gap),2)/loop_gap);hold on
    %loglog(KX_1d,-(sum(specGR_record(:,(loopindex-1)*loop_gap+1:loopindex*loop_gap),2)/loop_gap)/epsilon);hold on
    if(color_Num==2)
     figure(7);loglog(KX_1d,sum(specGG_record(:,(loopindex-1)*loop_gap+1:loopindex*loop_gap),2)/loop_gap);hold on
     figure(8);loglog(KX_1d,sum(specRR_record(:,(loopindex-1)*loop_gap+1:loopindex*loop_gap),2)/loop_gap);hold on
     figure(9);loglog(KX_1d,-(sum(specGR_record(:,(loopindex-1)*loop_gap+1:loopindex*loop_gap),2)/loop_gap));hold on
%      figure(7);plot(KX_1d,(sum(specGG_record(:,(loopindex-1)*loop_gap+1:loopindex*loop_gap),2)/loop_gap));hold on
%      figure(8);plot(KX_1d,(sum(specRR_record(:,(loopindex-1)*loop_gap+1:loopindex*loop_gap),2)/loop_gap));hold on
%      figure(9);plot(KX_1d,(sum(specGR_record(:,(loopindex-1)*loop_gap+1:loopindex*loop_gap),2)/loop_gap));hold on
     end
    legendseries = [legendseries, ["",num2str(loopindex*loop_gap)]];
end
if(0)  
 if(color_Num==1)
     if(velMode == 1)
         expected_specX = alpha*mean(c_noNoise(:,2:NX/2),1) + epsilon*avgr^2.*Ctwo(KX_1d*sigma)./(chi*KX_1d.^3 + epsilon*mean(c_noNoise(:,2:NX/2),1).*KX_1d.^2.*Cone(KX_1d*sigma));
     elseif(velMode == 2)
         expected_specX = alpha*mean(c_noNoise(:,2:NX/2),1) + epsilon*avgr^2.*Ctwo2(KX_1d*sigma,sigma)./(chi*KX_1d.^3 );
     end
 else
     if(velMode == 2)
     expected_specX = alpha*c_background*ones(size(KX_1d));
     expected_specGG = alpha*mean(c_noNoise(:,2:NX/2,1),1) + epsilon*avgr^2.*Ctwo2(KX_1d*sigma,sigma)./(chi*KX_1d.^3);
     expected_specRR = alpha*mean(c_noNoise(:,2:NX/2,2),1) + epsilon*avgr^2.*Ctwo2(KX_1d*sigma,sigma)./(chi*KX_1d.^3);
     expected_specGR = epsilon*avgr^2.*Ctwo2(KX_1d*sigma,sigma)./(chi*KX_1d.^3 );
     elseif(velMode == 1)
     expected_specX = alpha*c_background*ones(size(KX_1d));
     expected_specGG = alpha*mean(c_noNoise(:,2:NX/2,1),1) + epsilon*avgr^2.*Ctwo(KX_1d*sigma)./(chi*KX_1d.^3);
     expected_specRR = alpha*mean(c_noNoise(:,2:NX/2,2),1) + epsilon*avgr^2.*Ctwo(KX_1d*sigma)./(chi*KX_1d.^3);
     expected_specGR = epsilon*avgr^2.*Ctwo(KX_1d*sigma)./(chi*KX_1d.^3 );  
%      c_background_modified = c_background/2;
%      expected_specGG =  epsilon*c_background_modified^2.*Cone(KX_1d*sigma)./(chi*KX_1d);
%      expected_specRR =  epsilon*c_background_modified^2.*Cone(KX_1d*sigma)./(chi*KX_1d);
%      expected_specGR = epsilon*c_background_modified^2.*Cone(KX_1d*sigma)./(chi*KX_1d);

     end
 end
 
 figure(6);loglog(KX_1d,expected_specX);legend(legendseries(2:2:end),'expected');title('density spectrum');
 if(color_Num==2)
 figure(7);loglog(KX_1d,expected_specGG);legend(legendseries(2:2:end),'expected');title('green spectrum');hold off
 saveas(gcf,[filename,'/spec_GG.png'])
 figure(8);loglog(KX_1d,expected_specRR);legend(legendseries(2:2:end),'expected');title('red spectrum');hold off
 saveas(gcf,[filename,'/spec_RR.png'])
 figure(9);loglog(KX_1d,expected_specGR);legend(legendseries(2:2:end),'expected');title('red-green spectrum');hold off
 saveas(gcf,[filename,'/spec_GR.png'])
 end
end
if(0) %% plot data from Prof. Donev
% load('S_k_av.LFHD.1.dat')
% load('S_k_av.LFHD.2.dat')
% figure(6);
% loglog(KX_1d, S_k_av_LFHD_1(1:end-1,2)*epsilon);
% loglog(KX_1d, S_k_av_LFHD_2(1:end-1,2)*epsilon);
% loglog(KX_1d,epsilon*5.056862652e-6./KX_1d.^4);
% legend('PS data 1','PS data 2','FHD data 1','FHD data 2','theory')
load('Q2D_first.dat')
load('Q2D_second.dat')
figure(6);
loglog(Q2D_first(1:NX/2-1,1),Q2D_first(1:NX/2-1,2));
loglog(Q2D_second(1:NX/2-1,1),Q2D_second(1:NX/2-1,2));
%loglog(KX_1d,epsilon*5.056862652e-6./KX_1d.^4);
legend('PS data 1','PS data 2','FHD data 1','FHD data 2','theory')
end
xlabel('wave number'); ylabel('coefficient');hold off;
figure(6);
saveas(gcf,[filename,'/spec_X.png'])
end    





if (CYmean_gap && CYmean_plot)  %% movie for c(y) evolution

%%% recompute the initial

load([filename,'/CYmean.mat']);

CYmean = CYmean_record;
if(color_Num==2)
CYmean1 = CYmean_record_green;
CYmean2 = CYmean_record_red;
end

if(compare)
load([filename2,'/CYmean.mat']);
CYmean_comp = CYmean_record;
if(color_Num==2)
CYmean1_comp = CYmean_record_green;
CYmean2_comp = CYmean_record_red;
end
end

 for j = 1:length(CYmean_record(NY,:,1))
     figure(7);plot(CYmean(:,j));hold on;plot(sum(c_noNoise(:,1,:),3));
     if(color_Num==2)
     figure(8);plot(CYmean1(:,j));hold on;plot(c_noNoise(:,1,1));
     figure(9);plot(CYmean2(:,j));hold on;plot(c_noNoise(:,1,2));
     end
     if(compare)
        figure(7);plot(CYmean_comp(:,j)); 
        if(color_Num==2)
        figure(8);plot(CYmean1_comp(:,j)); 
        figure(9);plot(CYmean2_comp(:,j)); 
        end
     end
     figure(7);hold off;legend('full','initial','diffusion')
     title(['iter=',num2str((j-1)*CYmean_gap),' time=',num2str((j-1)*CYmean_gap*dt),' CFL = ',num2str(CFL)]);
     Movie(j)=getframe(gcf);
     if(color_Num==2)
     figure(8);hold off;legend('full','initial','diffusion')
     title(['iter=',num2str((j-1)*CYmean_gap),' time=',num2str((j-1)*CYmean_gap*dt),' CFL = ',num2str(CFL)]);
     Movie2(j)=getframe(gcf);
     figure(9);hold off;legend('full','initial','diffusion')
     title(['iter=',num2str((j-1)*CYmean_gap),' time=',num2str((j-1)*CYmean_gap*dt),' CFL = ',num2str(CFL)]);
     Movie3(j)=getframe(gcf);
     end
     %pause
 end
   
video = VideoWriter([filename,'/CY_whole.avi']); video.FrameRate=1;
open(video);writeVideo(video,Movie);close(video)
if(color_Num==2)
video = VideoWriter([filename,'/CY_G.avi']); video.FrameRate=1;
open(video);writeVideo(video,Movie2);close(video)
video = VideoWriter([filename,'/CY_R.avi']); video.FrameRate=1;
open(video);writeVideo(video,Movie3);close(video)
end
end








