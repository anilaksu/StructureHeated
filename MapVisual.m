% this routine is generated to visualize 2D solution
 
SimData1=importdata('2DNumerical.dat');
x=SimData1(:,1);
y=SimData1(:,2);
AmpFharx=SimData1(:,3);



Nx=20;
Ny=20;

%% let's reshape vectors to plot them

AmpFharx_r = reshape(AmpFharx,Nx,Ny);
x_r = reshape(x,Nx,Ny);
y_r = reshape(y,Nx,Ny);


%% let's generate the countor plot
figure(1)

contourf(x_r,y_r,AmpFharx_r,'LineColor','none')
colorbar
title('Induced Electric Field Within Solution Domain', 'FontSize', 10)
xlabel('x/L_k')
ylabel('y/L_k')
hcb=colorbar
set(gca,'fontsize',16)
title(hcb,'E/E_0')
pbaspect([2 1 1])
