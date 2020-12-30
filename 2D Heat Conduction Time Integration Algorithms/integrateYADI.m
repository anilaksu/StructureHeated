function [ Temp] = integrateYADI( Temp_x, Nx, Ny, dt, Q_yADI,S_yADI)
% this function performs ADI time integration in y direction

Temp = zeros((Ny-1)*(Nx-1),1); 
for i = 1:(Nx-1)
    for j = 1:(Ny-1)
         T_aux(j,1) = Temp_x(i+(j-1)*(Nx-1));  % vector mapping
    end
    %Temp((i-1)*(Ny-1)+1:i*(Ny-1),1) = S_yADI(:,:,i)*(Temp_x((i-1)*(Ny-1)+1:i*(Ny-1))+dt*Q_yADI(:,i));
    T_aux = S_yADI(:,:,i)*(T_aux+dt*Q_yADI(:,i));
    for j = 1:(Ny-1)
        Temp(i+(j-1)*(Nx-1))= T_aux(j,1);  % vector remapping
    end
end

end

