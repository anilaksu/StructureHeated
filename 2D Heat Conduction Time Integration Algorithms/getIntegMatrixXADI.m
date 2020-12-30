function [ S_xADI ] = getIntegMatrixXADI( K_xADI, Nx, Ny,dt )
% this function generates the solver matrix in x direction for ADI
S_xADI = zeros((Nx-1),(Nx-1),(Ny-1));

for i=1:(Ny-1)
    S_xADI(:,:,i)= inv(eye(Nx-1)-dt*K_xADI(:,:,i));
end


end

