% test_lawnmower_10step.m
%
% 10-iteration multi-step particle filter experiment: a sensor sweeps a
% lawnmower (boustrophedon) pattern across a STATIC plume, taking one
% measurement per iteration. Tests whether multiple, geometrically
% separated readings localise the source better than the single fixed
% reading of Report I Section 5.5 -- and exposes the static-source
% resampling cost (Section 5.3.3/5.4) as it actually accumulates over
% several steps rather than just one.
%
% CONFIGURATION (confirmed with Aryan):
%   Grid 60x60, dx=dy=1.0, source at (10, 30).
%   u, v, D, Q, dt carried over from the validated Day-10 configuration
%   (u=1.0, v=0.0, D=0.5, Q=1.0, dt=0.1) -- NOT independently re-specified
%   for this experiment, flagged here for visibility.
%
% LAWNMOWER PATH (one U-turn, 10 positions total):
%   Derived from the ACTUAL extracted level-set contour's bounding box
%   (not guessed from sigma -- see the comment at sensorPath below for why
%   that was wrong). Two columns at 30%/70% of the contour's x-extent,
%   5 y-positions each spanning the contour's y-extent plus a small margin:
%   column 1 sweeps up, U-turn, column 2 sweeps back down.
%
% RESAMPLING POLICY: existing ESS < N/2 trigger (Section 5.3.1), NOT
% forced every iteration. A ROUGHENING KERNEL (Gordon, Salmond & Smith
% 1993) is now applied immediately after every resampling event -- see
% ParticleFilterSIR.roughen(). This is a DELIBERATE methodological
% addition beyond what your report currently documents: Equation 46
% (theta_k = theta_{k-1}, zero process noise) no longer exactly describes
% what this script runs, since roughening injects small artificial noise
% specifically to counteract the permanent-diversity-loss structural gap
% Section 5.3.3/5.4 identified. Log this as a deliberate scope change in
% your diary/report, not as "the same static-source filter, just fixed."
%
% LEVEL SET OVERLAY: wired in via LevelSetExtractor, with the x/y
% index-convention swap fixed directly in that class (see LevelSetExtractor.m
% header comment for the full explanation).

clear; clc; close all;
rng(42);   % reproducibility, matches existing convention

% ---- output folder for .eps figures ----------------------------------------
% Matches the ./figs/ convention your LaTeX report's \stefig macro expects.
% Adjust the relative path if you run this script from somewhere other
% than tests/.
figDir = fullfile('..', 'figs', 'lawnmower');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

% ---- ground-truth forward model -------------------------------------------
Nx = 60; Ny = 60;
dx = 1.0; dy = 1.0;
dt = 0.1;
u  = 1.0; v = 0.0;
D  = 0.5;
Q  = 1.0;

src_i = 10;  src_j = 30;
% NOTE: at dx=dy=1, the flagged round(src_i/dx) vs (idx-1)*dx mismatch is
% numerically invisible (round(src_i/1) = src_i), so this is safe for THIS
% experiment specifically. Do not reuse src_i/src_j as raw indices if you
% ever change dx away from 1 -- re-check the sanity plot below first.

ad = AdvectionDiffusion2D(Nx, Ny, dx, dy, dt, u, v, D, src_i, src_j);

warmupSteps = 200;
for n = 1:warmupSteps
    ad = ad.step_loop(Q);   % value class -- must reassign every step
end

% ---- level set extractor ----------------------------------------------------
% The x/y axis swap previously here has been FIXED DIRECTLY in
% LevelSetExtractor.m (contourc's column/row convention vs. our i/j=x/y
% convention). No call-site workaround needed anymore -- .x and .y from
% extract() are correctly physical x and y. Re-run the Section 4 circular
% validation case after updating the class file, to confirm no regression
% (it should still pass -- that test is swap-symmetric either way, so it
% can't prove the fix, but it should still be unbroken).
lse = LevelSetExtractor(dx, dy);
Cth = 0.1 * max(ad.C(:));   % ADJUST: arbitrary 10% of peak; pick whatever
                            % boundary definition you actually want

% ---- sanity check BEFORE trusting anything downstream ----------------------
% If the bright spot is not sitting at (10, 30), stop here -- fix the
% index convention before looking at any of the 10 iterations below.
figure;
imagesc(ad.C'); set(gca, 'YDir', 'normal'); axis equal tight; colorbar;
hold on;
plot(src_i, src_j, 'r*', 'MarkerSize', 14, 'LineWidth', 2);
title('SANITY CHECK: settled field vs assumed source location (10,30)');
hold off;
exportgraphics(gcf, fullfile(figDir, 'sanity_check.eps'));

% ---- particle filter prior --------------------------------------------------
N = 1000;
pf = ParticleFilterSIR(N, 0, Nx, 0, Ny);

% ---- lawnmower sensor path, derived from the ACTUAL level-set footprint ----
% Previous version used fixed columns (x=20, x=40) sized off the Gaussian's
% sigma, which is the wrong quantity: sigma describes the underlying
% distribution's spread, not where a FIXED absolute threshold Cth actually
% crosses it. Since centreline concentration decays as 1/sqrt(x-x_s), the
% Cth-contour may not extend nearly as far downwind as sigma suggests --
% which is exactly why the old path drifted away from the boundary by the
% second column. Deriving the path from the real contour fixes this
% regardless of D, u, Q, or Cth, rather than re-guessing numbers by eye.
contour0 = lse.extract(ad.C, Cth);
if isempty(contour0)
    error('No level-set contour found at Cth=%.4f -- lower Cth or check the field.', Cth);
end
allX = horzcat(contour0.x);
allY = horzcat(contour0.y);

margin = 2;   % grid cells of slack outside the contour, so the sensor
              % travels just past the boundary rather than exactly on it
xLo = min(allX) - margin;  xHi = max(allX) + margin;
yLo = min(allY) - margin;  yHi = max(allY) + margin;

% Same snake pattern as before (up, U-turn, down), just confined to the
% contour's actual footprint instead of a guessed wide range.
xCol1 = xLo + 0.3*(xHi - xLo);
xCol2 = xLo + 0.7*(xHi - xLo);
yRows = linspace(yLo, yHi, 5);

sensorPath = [xCol1*ones(5,1), yRows'; ...
              xCol2*ones(5,1), fliplr(yRows)'];

fprintf('Tightened path bounding box: x in [%.1f, %.1f], y in [%.1f, %.1f]\n', ...
    xLo, xHi, yLo, yHi);

sigma_v = 0.01;
nIter = size(sensorPath, 1);

% ---- diagnostics storage ------------------------------------------------------
H_hist        = zeros(nIter, 1);
ESS_hist      = zeros(nIter, 1);
distinct_hist = zeros(nIter, 1);   % post-resample, PRE-roughening -- still
                                    % honestly shows the raw collapse (see note
                                    % at spread_hist below for why)
spread_hist   = zeros(nIter, 1);   % particle cloud std, AFTER roughening --
                                    % the meaningful diagnostic once roughening
                                    % is in play, since continuous jitter makes
                                    % distinct_hist trivially read N afterward
resampled     = false(nIter, 1);

fprintf('%4s %10s %10s %10s %10s\n', 'iter', 'ESS', 'H(nats)', 'distinct', 'resampled');

for k = 1:nIter
    sx = sensorPath(k, 1);
    sy = sensorPath(k, 2);
    sensor = PointSensor(sx, sy, sigma_v);
    y_k = sensor.measure(ad);

    pf.update(y_k, sensor, u, D, Q);   % handle class -- mutates in place

    % entropy/ESS MUST be read before resample (resample resets weights
    % to 1/N and would falsely report maximum entropy)
    H_hist(k)   = pf.computeEntropy();
    ESS_hist(k) = pf.computeESS();

    if ESS_hist(k) < N/2
        pf.resample();
        resampled(k) = true;

        % measured HERE, post-resample but PRE-roughening, so it still
        % honestly shows the raw collapse rather than being masked by jitter
        distinct_hist(k) = size(unique(pf.particles, 'rows'), 1);

        pf.roughen();   % K=0.2 default -- see ParticleFilterSIR.roughen()
    else
        distinct_hist(k) = size(unique(pf.particles, 'rows'), 1);
    end

    % spread, not distinct count, is the meaningful diagnostic once
    % roughening has fired even once -- continuous jitter makes every
    % particle numerically unique, so distinct_hist would trivially read
    % N forever after and tell you nothing about whether roughening
    % actually worked
    spread_hist(k) = std(pf.particles(:,1)) + std(pf.particles(:,2));

    fprintf('%4d %10.1f %10.4f %10d %10d\n', ...
        k, ESS_hist(k), H_hist(k), distinct_hist(k), resampled(k));

    % ---- plot this iteration (sqrt mode -- validated as the clearer encoding)
    plotParticleWeights(pf.particles, pf.weights, 'sqrt', ...
        struct('truePos', [src_i, src_j], 'sensorPos', [sx, sy], ...
               'title', sprintf('Iteration %d/%d: sensor (%d,%d)%s', ...
                   k, nIter, sx, sy, ternary(resampled(k), ' [RESAMPLED]', ''))));

    % overlay the path walked so far, so the lawnmower sweep is visible
    hold on;
    plot(sensorPath(1:k,1), sensorPath(1:k,2), 'k--', 'LineWidth', 1);
    hold off;

    % ---- level set overlay --------------------------------------------------
    contoursRaw = lse.extract(ad.C, Cth);
    hold on;
    for c = 1:numel(contoursRaw)
        plot(contoursRaw(c).x, contoursRaw(c).y, 'g-', 'LineWidth', 1.5);
    end
    hold off;

    exportgraphics(gcf, fullfile(figDir, sprintf('iteration_%d.eps', k)));
end

% ---- summary diagnostics across all 10 iterations -----------------------------
figure;
subplot(4,1,1);
plot(1:nIter, ESS_hist, 'o-'); yline(N/2, 'r--', 'N/2 trigger');
ylabel('ESS'); title('Diagnostics across 10 lawnmower iterations');

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

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end