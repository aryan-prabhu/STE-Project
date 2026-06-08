clear; close all;
dx = 1; dy = 1;
Lx = 40; Ly = 60;
Nx = round(Lx/dx); Ny = round(Ly/dy);
solver = AdvectionDiffusion2D(Nx,Ny,dx,dy,0.01,0.5,0.0,0.2,15,30);
solver.boundary_type = 'absorbing';
for k = 1:200000
    C_old = solver.C;
    solver = solver.step(10);
    residual = max(abs(solver.C(:)-C_old(:))) / max(abs(solver.C(:)));
    if residual < 1e-6, break, end
end
fprintf('stopped at step %d\n', k);
imagesc(solver.C'); axis xy; colorbar;
title('Absorbing — settled field');