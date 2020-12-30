function [ K_x, Q_x ] = getConductionX( T_1x, Q_2x, Nx,Ny, dx, dy )
% this function generates the conduction matrix in x direction with
% Dirichlet boundary condition on the left and Neumann Boundary condition
% on the right
% Nx: number of grid points
% dx: the step size in x direction 
% Q_x: the heat flow vector in x direction 

K_x1 = zeros ((Nx-1),(Nx-1));
Q_x1 = zeros ((Nx-1),1);

K_x = zeros ((Ny-1)*(Nx-1),(Ny-1)*(Nx-1));
Q_x = zeros ((Ny-1)*(Nx-1),1);

for i=1:(Nx-1)
   if (i == 1) % Dirichlet Boundary Condition 
        K_x1(i,i)= -2/(dx^2.);
        K_x1(i,(i+1))= 1/(dx^2.);
        Q_x1(i)= T_1x/(dx^2.);
   else if (i == (Nx-1)) % Neumann condition
           K_x1(i,(i-1)) = 1/(dx^2.);
           K_x1(i,i) = -1/(dx^2.);
           Q_x1(i)= Q_2x/dx;
       else
           K_x1(i,(i-1)) = 1/(dx^2.);
           K_x1(i,i) = -2/(dx^2.);
           K_x1(i,(i+1)) = 1/(dx^2.);
   end
end

for i = 1:(Ny-1)
   Y=eye(Nx-1)*(0.5+0.5*sqrt((i-1)*dy));
   K_x(((i-1)*(Nx-1)+1):i*(Nx-1),((i-1)*(Nx-1)+1):i*(Nx-1)) = Y*K_x1;
   Q_x(((i-1)*(Nx-1)+1):i*(Nx-1)) =  Y*Q_x1;
end

end

