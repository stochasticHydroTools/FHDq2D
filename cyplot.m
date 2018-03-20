
clear all;
filename = 'eps1vel1conv0noise0color_Num2VM3_kc0';
load([filename,'/parameters.mat']);
load([filename,'/CYmean.mat']);

section = floor(MaxIter/2);

CYmean_red_1_64 = mean(CYmean_record_red(:,2:section+1),2);
CYmean_red_2_64 = mean(CYmean_record_red(:,section+2:MaxIter+1),2);
CYmean_green_1_64 = mean(CYmean_record_green(:,2:section+1),2);
CYmean_green_2_64 = mean(CYmean_record_green(:,section+2:MaxIter+1),2);


 Y = LY/(NY-1)*(0:1:NY-1);
 
 CY_green_64_1 = [Y' CYmean_green_1_64];
 CY_green_64_2 = [Y' CYmean_green_2_64];
 CY_red_64_1 = [Y' CYmean_red_1_64];
 CY_red_64_2 = [Y' CYmean_red_2_64];
 
 save([filename,'/CY_green_64_eps1_kc0_1.txt'],'CY_green_64_1','-ascii')
 save([filename,'/CY_green_64_eps1_kc0_2.txt'],'CY_green_64_2','-ascii')
 
 save([filename,'/CY_red_64_eps1_kc0_1.txt'],'CY_red_64_1','-ascii')
 save([filename,'/CY_red_64_eps1_kc0_2.txt'],'CY_red_64_2','-ascii')
 
