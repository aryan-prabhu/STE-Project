clear; close all;
Nx=40; Ny=60; dx=1; dy=1; dt=0.01;
u=0.5; v=0.0; D=0.2; Q=10;
src_x=15; src_y=30;

solver = AdvectionDiffusion2D(Nx,Ny,dx,dy,dt,u,v,D,src_x,src_y);

tol       = 1e-6;
max_steps = 100000;
residuals = zeros(max_steps,1);

for k = 1:max_steps
    C_old  = solver.C;
    solver = solver.step(Q);
    C_new  = solver.C;

    residual = max(abs(C_new(:) - C_old(:))) / max(abs(C_new(:)));
    residuals(k) = residual;

    if residual < tol
        fprintf('Steady state reached at step %d (residual %.2e)\n', k, residual);
        residuals = residuals(1:k);   % trim unused tail
        break
    end
end
% --- compare settled solver field against analytical Gaussian plume ---
src_i = solver.src_i;  src_j = solver.src_j;
C_analytic = gaussianPlume(Nx,Ny,dx,dy,u,D,Q,src_i,src_j);
C_num = solver.C;

% (the old unmasked L2_error line — fine to leave it in OR comment it out,
%  it's harmless either way)
L2_error = norm(C_num(:) - C_analytic(:)) / norm(C_analytic(:));
fprintf('Relative L2 error: %.4f  (%.2f%%)\n', L2_error, 100*L2_error);

% --- masked error: exclude a neighbourhood of the singular source ---
exclude_radius = 5;
[I, J] = ndgrid(1:Nx, 1:Ny);
dist_to_src = sqrt((I - src_i).^2 + (J - src_j).^2);
mask = dist_to_src > exclude_radius;

err_field = (C_num - C_analytic) .* mask;
ref_field =  C_analytic          .* mask;

L2_error_masked = norm(err_field(:)) / norm(ref_field(:));
fprintf('Masked relative L2 error (r>%d): %.4f  (%.2f%%)\n', ...
        exclude_radius, L2_error_masked, 100*L2_error_masked);



figure;
subplot(1,3,1); imagesc(C_num');      axis xy; colorbar; title('Numerical');
subplot(1,3,2); imagesc(C_analytic'); axis xy; colorbar; title('Analytical');
subplot(1,3,3); imagesc((C_num-C_analytic)'); axis xy; colorbar; title('Difference');



figure; semilogy(residuals); grid on;
xlabel('step'); ylabel('relative residual');
title('Convergence to steady state');