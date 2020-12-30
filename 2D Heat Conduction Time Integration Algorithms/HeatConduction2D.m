% This script is developed to solve 2-D Conduction Problem with variable
% coefficients


clear all 
format long

Nx= 6;          % Number of grid points in x direction
Ny= 3;          % Number of grid points in y direction
N_time = 400;   % Number of time steps
dx = 1;         % the step size in x direction
dy = 1;         % the step size in y direction
dt = 1;         % the time step size 


%% Boundary Conditions
Tx = 623;       % Dirichlet condition boundary temperature in x direction
Ty = 473;       % Dirichlet condition boundary temperature in y direction
Qx = 3;         % Neumann condition boundary heat flow in x direction
Qy = -4;        % Neumann condition boundary heat flow in y direction

Temp_Crank=zeros((Nx-1)*(Ny-1),N_time);   % Temperature array for Crank-Nicholson
Temp_ADI=zeros((Nx-1)*(Ny-1),N_time);   % Temperature array for Crank-Nicholson

Temp_Crank(:,1)= 300;  % Initial Condition 
Temp_ADI(:,1)= 300;  % Initial Condition 

[ K_x, Q_x ] = getConductionX( Tx, Qx, Nx, Ny, dx, dy );         % Classical Conduction Matrix in x direction
[ K_xADI, Q_xADI ] = getConductionXADI( Tx, Qx, Nx,Ny, dx, dy ); % Alternating Direction Implicit in x direction
[ S_xADI ] = getIntegMatrixXADI( K_xADI, Nx, Ny,dt );            % Solver matrices in x direction for ADI

[ K_y, Q_y ] = getConductionY( Ty, Qy, Nx, Ny, dx, dy );         % Classical Conduction Matrix in y direction
[ K_yADI, Q_yADI ] = getConductionYADI( Ty, Qy, Nx,Ny, dx, dy ); % Alternating Direction Implicit in y direction
[ S_yADI ] = getIntegMatrixYADI( K_yADI, Nx, Ny,dt );            % Solver matrices in y direction for ADI

K = K_x + K_y; %  total conduction matrix
Q = Q_x + Q_y; %  Total Heat flow vector

% time integration 
tic
M_crank = eye((Nx-1)*(Ny-1))-0.5*dt*K;  % Matrix multiplied with T_(n+1)
M_inv = inv(M_crank);
for i=2:N_time
   % Crank-Nicholson Time Integration
   Temp_Crank(:,i) =  M_inv*(0.5*dt*K+eye((Nx-1)*(Ny-1)))*Temp_Crank(:,i-1)+M_inv*dt*Q;
   % ADI time Integration in x direction
   [ Temp_ADI(:,i) ] = integrateXADI( Temp_ADI(:,i-1), Nx, Ny, dt, Q_xADI,S_xADI);
   % ADI time Integration in y direction
   [ Temp_ADI(:,i)] = integrateYADI( Temp_ADI(:,i), Nx, Ny, dt, Q_yADI,S_yADI);
end
toc

% Boundary Temperatures added to Temperature distribution found by Crank
% Nicholson
[ Temp ] = addTempBoundaries( Temp_Crank(:,N_time), Nx, Ny, Tx, Ty, Qx, Qy, dx, dy  );

x = linspace(0,6,7);
y = linspace(0,3,4);

%Temp_final = reshape(Temp_aux(:,N_time),[5, 2]);
TempCrank_final = reshape(Temp,[7, 4]);
set(0,'DefaultAxesFontSize',13)
figure(1)
contourf(x,y,transpose(TempCrank_final)), shading flat
colormap('jet');
colorbar;
% %caxis([-1, 1]);
title('Temperature Distribution T(x,y) with Crank Nicholson') % title
ylabel('y')
% % label for y axis
xlabel('x')
% pbaspect([1 1 0.2])

% Boundary Temperatures added to Temperature distribution found by Crank
% Nicholson
[ Temp ] = addTempBoundaries( Temp_ADI(:,N_time), Nx, Ny, Tx, Ty, Qx, Qy, dx, dy  );

TempADI_final = reshape(Temp,[7, 4]);
set(0,'DefaultAxesFontSize',13)
figure(2)
contourf(x,y,transpose(TempADI_final)), shading flat
colormap('jet');
colorbar;
% %caxis([-1, 1]);
title('Temperature Distribution T(x,y) with ADI') % title
ylabel('y')
% % label for y axis
xlabel('x')
% pbaspect([1 1 0.2])
