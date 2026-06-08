clear classes; close all;

Nx=40; Ny=60; dx=1; dy=1; dt=0.01;
u=0.5; v=0.0; D=0.2; src_x=15; src_y=30;

solverA = AdvectionDiffusion2D(Nx,Ny,dx,dy,dt,u,v,D,src_x,src_y);
solverB = AdvectionDiffusion2D(Nx,Ny,dx,dy,dt,u,v,D,src_x,src_y);

for k = 1:200
    solverA = solverA.step_loop(10);   % old
    solverB = solverB.step(10);        % new vectorised
end

max_abs_diff = max(abs(solverA.C(:) - solverB.C(:)))