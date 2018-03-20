function [cMean_record, spec2dmean] = real_spec_2dmean(c, chat, iter, real_spec_2dMean_gap, NX, NY, color_Num, cMean_record, spec2dmean)
% computes average of mean of S(kx,ky) over kx and ky

%%%% mean for each color
for j = 1:color_Num
     cMean_record(iter/real_spec_2dMean_gap+1,j) = cMean_record(iter/real_spec_2dMean_gap+1,j) + mean(mean(c(:,:,j)));
     
     spec2d = chat(:,:,j).*conj(chat(:,:,j));
     spec2d_nzwave = spec2d(2:NY,2:NX); 
     spec2dmean(iter/real_spec_2dMean_gap+1,j) = spec2dmean(iter/real_spec_2dMean_gap+1,j) + mean(spec2d_nzwave(:));   
end
%%% mean for whole
cMean_record(iter/real_spec_2dMean_gap+1,color_Num+1) = cMean_record(iter/real_spec_2dMean_gap+1,color_Num+1) + mean(c(:));
     
chat_whole = sum(chat,3);
spec2d = chat_whole.*conj(chat_whole);
spec2d_nzwave = spec2d(2:NY,2:NX); 
spec2dmean(iter/real_spec_2dMean_gap+1,color_Num+1) = spec2dmean(iter/real_spec_2dMean_gap+1,color_Num+1) + mean(spec2d_nzwave(:)); 
end

