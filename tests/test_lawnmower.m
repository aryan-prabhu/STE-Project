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
%   Column 1 at x=20 (10 units downwind of source, sigma~3.2 there):
%       y = 10, 20, 30, 40, 50    (sweeping up)
%   Column 2 at x=40 (30 units downwind, sigma~5.5 there):
%       y = 50, 40, 30, 20, 10    (sweeping back down)
%   Column spacing was chosen using sigma^2 = 2*D*x/u so both columns
%   actually cross the plume core and both edges, rather than guessing.
%
% RESAMPLING POLICY: existing ESS < N/2 trigger (Section 5.3.1), NOT
% forced every iteration. Because the source is a static parameter
% (Equation 46, theta_k = theta_{k-1}), there is no prediction step to
% perturb particles apart between resamples -- any diversity lost at a
% resampling event is PERMANENT (Section 5.3.3/5.4, your own documented
% structural gap). Watch the "distinct" column printed below: if it
% collapses hard early and stays flat, that is the gap your report
% already predicted, not a bug in this script.
%
% LEVEL SET OVERLAY: *** NOT YET WIRED IN -- see marked block below ***
% LevelSetExtractor's exact constructor/method signature is not in my
% notes (only that it wraps contourc and does index<->physical
% conversion). Given the ALREADY-FLAGGED index-convention mismatch
% between AdvectionDiffusion2D and LevelSetExtractor, guessing a third
% convention here risks silently misplacing the contour. Paste
% LevelSetExtractor.m and this gets wired in correctly.

clear; clc; close all;
rng(42);   % reproducibility, matches existing convention

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

% ---- particle filter prior --------------------------------------------------
N = 1000;
pf = ParticleFilterSIR(N, 0, Nx, 0, Ny);

% ---- lawnmower sensor path (10 positions, one U-turn) ------------------------
sensorPath = [20, 10; 20, 20; 20, 30; 20, 40; 20, 50; ...
              40, 50; 40, 40; 40, 30; 40, 20; 40, 10];
sigma_v = 0.01;
nIter = size(sensorPath, 1);

% ---- diagnostics storage ------------------------------------------------------
H_hist        = zeros(nIter, 1);
ESS_hist      = zeros(nIter, 1);
distinct_hist = zeros(nIter, 1);
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
    end

    % distinct count read AFTER the resample decision, so it reflects
    % this iteration's actual outcome (mirrors Table 3's "post-resample"
    % framing)
    distinct_hist(k) = size(unique(pf.particles, 'rows'), 1);

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
end

% ---- summary diagnostics across all 10 iterations -----------------------------
figure;
subplot(3,1,1);
plot(1:nIter, ESS_hist, 'o-'); yline(N/2, 'r--', 'N/2 trigger');
ylabel('ESS'); title('Diagnostics across 10 lawnmower iterations');

subplot(3,1,2);
plot(1:nIter, H_hist, 'o-');
ylabel('Entropy (nats)');

subplot(3,1,3);
plot(1:nIter, distinct_hist, 'o-');
ylabel('Distinct particles'); xlabel('iteration');

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end