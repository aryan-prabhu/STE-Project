clear; clc; close all;
rng(42);

figDir = fullfile('..', 'figs', 'lawnmower');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

Nx = 60; Ny = 60;
dx = 1.0; dy = 1.0;
dt = 0.1;
u  = 1.0; v = 0.0;
D  = 0.5;
Q  = 1.0;

src_i = 10;  src_j = 30;

ad = AdvectionDiffusion2D(Nx, Ny, dx, dy, dt, u, v, D, src_i, src_j);

warmupSteps = 200;
for n = 1:warmupSteps
    ad = ad.step_loop(Q);
end

lse = LevelSetExtractor(dx, dy);
Cth = 0.1 * max(ad.C(:));

figure;
imagesc(ad.C'); set(gca, 'YDir', 'normal'); axis equal tight; colorbar;
hold on;
plot(src_i, src_j, 'r*', 'MarkerSize', 14, 'LineWidth', 2);
hold off;
exportgraphics(gcf, fullfile(figDir, 'sanity_check.eps'));

N = 1000;
theta_Q = 0.25;
pf = ParticleFilterSIR(N, 0, Nx, 0, Ny, theta_Q);

contour0 = lse.extract(ad.C, Cth);
if isempty(contour0)
    error('No level-set contour found at Cth=%.4f -- lower Cth or check the field.', Cth);
end
allX = horzcat(contour0.x);
allY = horzcat(contour0.y);

margin = 2;
xLo = min(allX) - margin;  xHi = max(allX) + margin;

xUpwind = src_i - margin;

yDownLo = src_j - 13;  yDownHi = src_j + 13;
yUpLo   = 0.5*Ny - 0.3*Ny;  yUpHi   = 0.5*Ny + 0.3*Ny;

yUpwind = linspace(yUpLo, yUpHi, 2)';
yRows6  = linspace(yDownLo, yDownHi, 6)';

xCols = xLo + [0.25, 0.50, 0.75] * (xHi - xLo) + 5;

sensorPath = [xUpwind*ones(2,1), yUpwind; ...
              xCols(1)*ones(6,1), yRows6; ...
              xCols(2)*ones(6,1), flipud(yRows6); ...
              xCols(3)*ones(6,1), yRows6];

fprintf('Path: upwind (x=%.1f), then 25%%->50%%->75%% downwind columns (x=%.1f,%.1f,%.1f), 20 steps total\n', ...
    xUpwind, xCols(1), xCols(2), xCols(3));

sigma_v = 0.01;
nIter = size(sensorPath, 1);

H_hist        = zeros(nIter, 1);
ESS_hist      = zeros(nIter, 1);
distinct_hist = zeros(nIter, 1);
spread_hist   = zeros(nIter, 1);
resampled     = false(nIter, 1);

fprintf('%4s %10s %10s %10s %10s\n', 'iter', 'ESS', 'H(nats)', 'distinct', 'resampled');

for k = 1:nIter
    sx = sensorPath(k, 1);
    sy = sensorPath(k, 2);
    sensor = PointSensor(sx, sy, sigma_v);
    y_k = sensor.measure(ad);

    pf.update(y_k, sensor, u, D);

    H_hist(k)   = pf.computeEntropy();
    ESS_hist(k) = pf.computeESS();

    if ESS_hist(k) < N/2
        fprintf('Pre-resample: %d particles with weight in top 1%%, their x_s range [%.1f, %.1f]\n', ...
        sum(pf.weights > prctile(pf.weights,99)), ...
        min(pf.particles(pf.weights > prctile(pf.weights,99),1)), ...
        max(pf.particles(pf.weights > prctile(pf.weights,99),1)));
        pf.resample();
        resampled(k) = true;

        distinct_hist(k) = size(unique(pf.particles, 'rows'), 1);

        pf.roughen();
    else
        distinct_hist(k) = size(unique(pf.particles, 'rows'), 1);
    end

    spread_hist(k) = std(pf.particles(:,1)) + std(pf.particles(:,2));

    fprintf('%4d %10.1f %10.4f %10d %10d\n', ...
        k, ESS_hist(k), H_hist(k), distinct_hist(k), resampled(k));

    plotParticleWeights(pf.particles, pf.weights, 'sqrt', ...
        struct('truePos', [src_i, src_j], 'sensorPos', [sx, sy]));

    hold on;
    plot(sensorPath(1:k,1), sensorPath(1:k,2), 'k--', 'LineWidth', 1);
    hold off;

    contoursRaw = lse.extract(ad.C, Cth);
    hold on;
    for c = 1:numel(contoursRaw)
        plot(contoursRaw(c).x, contoursRaw(c).y, 'g-', 'LineWidth', 1.5);
    end
    hold off;

    exportgraphics(gcf, fullfile(figDir, sprintf('iteration_%d.eps', k)));
end

figure;
subplot(4,1,1);
plot(1:nIter, ESS_hist, 'o-'); yline(N/2, 'r--', 'N/2 trigger');
ylabel('ESS');

subplot(4,1,2);
plot(1:nIter, H_hist, 'o-');
ylabel('Entropy (nats)');

subplot(4,1,3);
plot(1:nIter, distinct_hist, 'o-');
ylabel('Distinct (pre-roughen)');

subplot(4,1,4);
plot(1:nIter, spread_hist, 'o-');
ylabel('Cloud spread'); xlabel('iteration');

exportgraphics(gcf, fullfile(figDir, 'diagnostics.eps'));