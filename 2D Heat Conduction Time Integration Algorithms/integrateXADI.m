function [ Temp_x ] = integrateXADI( Temp_0, Nx, Ny, dt, Q_xADI,S_xADI)
% this function performs ADI time integration in x direction

Temp_x = zeros((Ny-1)*(Nx-1),1); 
for i = 1:(Ny-1)
    Temp_x((i-1)*(Nx-1)+1:i*(Nx-1),1) = S_xADI(:,:,i)*(Temp_0((i-1)*(Nx-1)+1:i*(Nx-1))+dt*Q_xADI(:,i));
end

end

