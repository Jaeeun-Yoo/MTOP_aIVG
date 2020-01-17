function plot_density(x,nelx,nely,XL,YL,loop,Mnd2)


% PLOT DENSITIES  
x_axis=(1:1:nelx)*XL/nelx;
y_axis=(1:1:nely)*YL/nely;
figure(1);
colormap(gray); 
imagesc(x_axis,y_axis,1-x); 
axis equal; 
axis tight; 
axis on;
%colorbar
%drawnow;

title(['Density plot (',num2str(loop),' iter. Gray level indicator = ',num2str(Mnd2,'%0.2f'),'% )'])
xlabel('x (meter)')
ylabel('y (meter)')

set(gcf,'position',[100 100 800 600]);
