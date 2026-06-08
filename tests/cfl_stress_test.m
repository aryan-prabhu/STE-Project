clear classes; close all;

Nx=40; Ny=60; dx=1; dy=1; dt=0.01;
v=0.0; D=0.2; src_x=15; src_y=30;

u_values = [0.5, 2, 5, 20, 100];   % climb until it breaks

for u = u_values
    solver = AdvectionDiffusion2D(Nx,Ny,dx,dy,dt,u,v,D,src_x,src_y);
    for k = 1:200
        solver = solver.step(10);
    end
    courant = u*dt/dx;
    peak    = max(solver.C(:));
    fprintf('u = %6.1f | Courant = %6.3f | peak C = %g\n', u, courant, peak);
end

fprintf('\n--- Diffusion sweep ---\n');
u = 0.5;
D_values = [0.2, 5, 20, 25, 40];

for D = D_values
    solver = AdvectionDiffusion2D(Nx,Ny,dx,dy,dt,u,v,D,src_x,src_y);
    for k = 1:200
        solver = solver.step(10);
    end
    diff_number = D*dt*(1/dx^2 + 1/dy^2);
    peak        = max(solver.C(:));
    fprintf('D = %5.1f | diff number = %6.3f | peak C = %g\n', D, diff_number, peak);
end