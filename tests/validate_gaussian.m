clear classes; close all;
Nx=40; Ny=60; dx=1; dy=1; u=0.5; D=0.2; Q=10;
src_i=15; src_j=30;

C_analytic = gaussianPlume(Nx,Ny,dx,dy,u,D,Q,src_i,src_j);

figure; imagesc(C_analytic'); axis xy; colorbar;
title('Analytical Gaussian plume');