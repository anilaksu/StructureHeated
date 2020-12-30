function [ K_y, Q_y ] = getConductionY( T_2y, Q_1y, Nx,Ny, dx, dy )
% this function generates the conduction matrix in y direction with
% Dirichlet boundary condition on the left and Neumann Boundary condition
% on the right
% Nx: number of grid points
% dx: the step size in x direction 
% Q_x: the heat flow vector in x direction 

K_y1 = zeros ((Ny-1),(Ny-1));
Q_y1 = zeros ((Ny-1),1);

K_y = zeros ((Ny-1)*(Nx-1),(Ny-1)*(Nx-1));
Q_y = zeros ((Ny-1)*(Nx-1),1);

for i=1:(Ny-1)
   if (i == (Ny-1)) % Dirichlet Boundary Condition 
        K_y1(i,i)= -2/(dy^2.);
        K_y1(i,(i-1))= 1/(dy^2.);
        Q_y1(i)= T_2y/(dy^2.);
   else if (i == 1) % Neumann Condition
           K_y1(i,i) = -1/(dy^2.);
           K_y1(i,(i+1)) = 1/(dy^2.);
           Q_y1(i)= -Q_1y/dy;
       else
           K_y1(i,(i-1)) = 1/(dy^2.);
           K_y1(i,i) = -2/(dy^2.);
           K_y1(i,(i+1)) = 1/(dy^2.);
   end
end

% for i = 1:(Nx-1)
%    X=eye(Ny-1)*(0.5+0.5*sqrt(i*dx));   % Variable coefficient
%    K_y(((i-1)*(Ny-1)+1):i*(Ny-1),((i-1)*(Ny-1)+1):i*(Ny-1)) = X*K_y1;
%    Q_y(((i-1)*(Ny-1)+1):i*(Ny-1)) = Q_y1;
% end

for i = 1:(Nx-1)
    X=eye(Ny-1)*(0.5+0.5*sqrt(i*dx));   % Variable coefficient
    K_dummy = X*K_y1 ;                  % Dummy conduction matrix to be update at each time step
    Q_dummy = X*Q_y1 ;                  % Dummy conduction matrix to be update at each time step
    for j = 1:(Ny-1)
        for k = 1:(Ny-1)
            K_y((j-1)*(Nx-1)+i,i+(k-1)*(Nx-1)) = K_dummy(j,k);
        end
        Q_y((j-1)*(Nx-1)+i) = Q_dummy(j);
    end 
end

end
