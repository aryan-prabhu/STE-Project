% test_plotParticleWeights.m
%
% Validates plotParticleWeights.m against the single-step case of
% Report I, Section 5.5: true source (3,4), sensor (6,4), N = 1000,
% uniform prior over a 10 x 8 domain.
%
% Run this BEFORE reusing plotParticleWeights in the 10-iteration
% lawnmower experiment. The qualitative answer here is already known
% (the mirror-symmetric ridge of Figure 7) -- this script exists to
% check that 'sqrt' size-encoding mode reproduces the SAME qualitative
% posterior shape as the existing 'color' mode, before trusting
% size-encoding on a result you don't already know the answer to.
% (A 'log' mode was tried and removed -- it blurred the ridge boundary
% by giving the near-dead tail too much of the size range. See
% plotParticleWeights.m header for the full explanation.)
%
% PLACEHOLDER VALUES TO VERIFY: sigma_v, Q, u, D below are NOT stated in
% the report text for this specific run -- only y_k = 0.230329 is (Table
% 3). The other four are carried over from the Day-10 experiment as a
% best guess. Check them against your actual Section-5.5 test script
% before trusting this output; a mismatch here would change the posterior
% shape and could look like a bug in the new plotting code when it isn't.

clear; clc; close all;
rng(42);   % reproducibility, matches existing convention

% ---- domain / true source / sensor (matches Section 5.5) --------------
x_min = 0;  x_max = 10;
y_min = 0;  y_max = 8;
truePos   = [3, 4];
sensorPos = [6, 4];

% ---- VERIFY THESE FOUR AGAINST YOUR ACTUAL SECTION 5.5 SCRIPT ----------
sigma_v = 0.01;
Q = 1.0;
u = 1.0;
D = 0.5;
% -------------------------------------------------------------------------

N = 1000;

% ---- particle filter prior ---------------------------------------------
pf = ParticleFilterSIR(N, x_min, x_max, y_min, y_max);

% ---- sensor and the one noisy reading -----------------------------------
sensor = PointSensor(sensorPos(1), sensorPos(2), sigma_v);
y_k = 0.230329;   % exact reading reported in Table 3 -- reuse it so this
                  % run is directly comparable to Figure 7 / Table 3

% ---- prior plot (all particles equal weight) -----------------------------
w_prior = ones(N, 1) / N;
modes = {'color', 'sqrt'};
for m = 1:numel(modes)
    plotParticleWeights(pf.particles, w_prior, modes{m}, ...
        struct('truePos', truePos, 'sensorPos', sensorPos, ...
               'title', sprintf('Prior (%s)', modes{m})));
end

% ---- one Bayesian update --------------------------------------------------
% ParticleFilterSIR is a handle class: this call mutates pf in place,
% no reassignment needed.
pf.update(y_k, sensor, u, D, Q);

% ---- posterior plot, before resampling ------------------------------------
% (entropy/ESS are not computed here -- this script is for the visual
% check only, not the diagnostic table)
for m = 1:numel(modes)
    plotParticleWeights(pf.particles, pf.weights, modes{m}, ...
        struct('truePos', truePos, 'sensorPos', sensorPos, ...
               'title', sprintf('Posterior (%s)', modes{m})));
end