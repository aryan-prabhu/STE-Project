clear classes
close all

% Build the solver
solver = AdvectionDiffusion2D(60, 60, 1.0, 1.0, 0.05, 6.0, 0.0, 2.0, 15.0, 30.0);

% Run 200 steps with constant source strength Q = 10
for k = 1:200
    solver = solver.step(10.0);
end

% Visualise the final concentration field
imagesc(solver.C')
axis xy
colorbar
title('Concentration field after 200 steps')