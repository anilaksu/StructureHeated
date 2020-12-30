function [ Temp ] = addTempBoundaries( T_aux, Nx, Ny, Tx, Ty, Qx, Qy, dx, dy  )
% this function adds boundaries to the computed temperature field

Temp_y = zeros((Ny+1)*(Nx-1),1); % let's first add boundary points in y direction
Temp_y(1:(Nx-1),1) =  T_aux(1:(Nx-1),1)-dy*Qy; % Neumann boundaries 
Temp_y((Nx-1)+1:(Ny)*(Nx-1),1) =  T_aux;       % Domain points
Temp_y(Ny*(Nx-1)+1:end,1) =  Ty;               % Dirichlet boundaries 

Temp = zeros((Ny+1)*(Nx+1),1); % let's add boundary points in x direction

for i=1:(Ny+1)
    Temp((i-1)*(Nx+1)+1,1)= Tx;                                          % Dirichlet boundaries 
    Temp((i-1)*(Nx+1)+2:i*(Nx+1)-1,1)= Temp_y((i-1)*(Nx-1)+1:i*(Nx-1),1) % Domain points
    Temp(i*(Nx+1),1)= Temp_y(i*(Nx-1),1)+ dx*Qx;                       % Neumann boundaries
end
end

