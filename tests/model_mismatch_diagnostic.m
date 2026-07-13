% test_model_mismatch_centerline.m
%
% Re-run with warmupSteps increased from 200 to 600 to test the
% steady-state hypothesis: at 200 steps (t=20), the advective front
% (u=1) had only reached x~30 by warm-up's end, meaning the far-field
% points in the previous sweep were comparing plumeConc (assumes full
% steady state) against a still-developing PDE transient -- not a
% structural model defect. At 600 steps (t=60), the front should have
% crossed the entire domain (needs t >= 58 for the sweep's farthest
% point), so the field should be genuinely settled everywhere sampled.

clear; clc; close all;
rng(42);

figDir = fullfile('..', 'figs', 'model_mismatch');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

Nx = 60; Ny = 60;
dx = 1.0; dy = 1.0;
dt = 0.1;
u  = 1.0; v = 0.0;
D  = 0.5;
Q  = 1.0;
src_i = 10; src_j = 30;

ad = AdvectionDiffusion2D(Nx, Ny, dx, dy, dt, u, v, D, src_i, src_j);

warmupSteps = 600;   % increased from 200 -- see header
for n = 1:warmupSteps
    ad = ad.step_loop(Q);
end

% ---- centerline sweep: y FIXED at src_j, x varies (unchanged) --------------
nPts  = 18;
xVals = linspace(src_i + 1, (Nx-2)*dx, nPts)';
yVals = src_j * ones(nPts, 1);

c_analytic = zeros(nPts, 1);
c_true     = zeros(nPts, 1);
ratio      = nan(nPts, 1);
pctError   = nan(nPts, 1);
downwindDist = xVals - src_i;

sensor0 = PointSensor(0, 0, 0);

for k = 1:nPts
    sx = xVals(k);
    sy = yVals(k);

    c_analytic(k) = plumeConc(src_i, src_j, sx, sy, u, D, Q);

    sensor0.x = sx;
    sensor0.y = sy;
    c_true(k) = sensor0.sampleField(ad);

    if abs(c_true(k)) > 1e-12
        ratio(k)    = c_analytic(k) / c_true(k);
        pctError(k) = 100 * (c_analytic(k) - c_true(k)) / c_true(k);
    end
end

fprintf('%6s %8s %12s %12s %10s %10s\n', ...
    'pt', 'dist', 'PDE', 'analytic', 'ratio', 'pctErr');
for k = 1:nPts
    fprintf('%6d %8.1f %12.4e %12.4e %10.4f %10.1f\n', ...
        k, downwindDist(k), c_true(k), c_analytic(k), ratio(k), pctError(k));
end

figure;
plot(downwindDist, pctError, 'o-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
xlabel('downwind distance from source (x - x_s), centerline');
ylabel('percent error: (analytic - PDE) / PDE \times 100');
yline(0, 'k--');
grid on;
exportgraphics(gcf, fullfile(figDir, 'pct_error_centerline_warmup600.eps'));
