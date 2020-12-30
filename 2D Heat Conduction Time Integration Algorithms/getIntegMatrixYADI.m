function [ S_yADI ] = getIntegMatrixYADI( K_yADI, Nx, Ny,dt )
% this function generates the solver matrix in y direction for ADI
S_yADI = zeros((Ny-1),(Ny-1),(Nx-1));

for i=1:(Nx-1)
    S_yADI(:,:,i)= inv(eye(Ny-1)-dt*K_yADI(:,:,i));
end


end
