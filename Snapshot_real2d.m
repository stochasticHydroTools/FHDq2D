function Movie = Snapshot_real2d(c, c_whole, clims, iter, dt, real2d_gap, color_Num, Movie)

for j = 1:color_Num
        colormap jet;
        figure(j);clf
        imagesc(c(:,:,j),clims);shading flat;colorbar;title(['iter=',num2str(iter),' time=',num2str(iter*dt),' dt = ',num2str(dt)]);
        Movie(iter/real2d_gap+1,j)=getframe(gcf);
end
colormap jet;
figure(color_Num+1);clf
imagesc(c_whole,clims);shading flat;colorbar;title(['iter=',num2str(iter),' time=',num2str(iter*dt),' dt = ',num2str(dt)]);
Movie(iter/real2d_gap+1,color_Num+1)=getframe(gcf);

end

