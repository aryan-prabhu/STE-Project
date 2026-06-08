%% test_update.m  --  Day 8: single measurement update, noise-free
% Goal: fire ParticleFilterSIR.update once and confirm the weights
% change sensibly. y_k is NOISE-FREE for this first run.

clear; clc;
%% --- physics constants (shared by truth and filter) ---
u = 1.0;        % wind speed in +x
D = 0.5;        % isotropic diffusion coefficient
Q = 1.0;        % source strength (known constant, Day 8)

%% --- the TRUTH (known to us, not to the filter) ---
x_s_true = 3.0;     % true source x
y_s_true = 4.0;     % true source y

%% --- the sensor (downwind of the true source so it smells the plume) ---
sensor = PointSensor(6.0, 4.0, 0.05);   % x, y, sigma

%% --- generate the noise-free reading y_k ---
% plumeConc gives the clean concentration; no sigma*randn added.
y_k = plumeConc(x_s_true, y_s_true, sensor.x, sensor.y, u, D, Q);
fprintf('Noise-free reading y_k = %.6f\n', y_k);

%% --- build the filter (Day 7 prior: uniform cloud over the domain) ---
N = 1000;
pf = ParticleFilterSIR(N, 0, 10, 0, 8);   % N, x_min, x_max, y_min, y_max
H_prior = pf.computeEntropy();   % prior entropy, before any update

w_before = pf.weights;                    % snapshot: should be uniform
%% --- PLOT: the prior cloud (before any update) ---
figure;
scatter(pf.particles(:,1), pf.particles(:,2), 25, w_before, 'filled');
colorbar; colormap(parula); hold on;
plot(x_s_true, y_s_true, 'rp', 'MarkerSize', 18, 'MarkerFaceColor', 'r');
plot(sensor.x,  sensor.y,  'ks', 'MarkerSize', 12, 'MarkerFaceColor', 'k');
xlabel('x'); ylabel('y');
title('Prior: uniform particle cloud (before update)');
axis([0 10 0 8]); axis equal; hold off;
caxis([0 7e-3]);   % shared colour scale (see note)
exportgraphics(gcf, 'prior_cloud.png', 'Resolution', 200);
%% --- fire the single update (Eq 4) ---
pf.update(y_k, sensor, u, D, Q);
w_after = pf.weights;
H_post  = pf.computeEntropy();   % capture BEFORE any resample

%% --- CHECKS ---
fprintf('--- checks ---\n');
fprintf('weights sum to 1     : %.10f\n', sum(w_after));
fprintf('all weights >= 0     : %d\n', all(w_after >= 0));
fprintf('variance before      : %.3e\n', var(w_before));
fprintf('variance after       : %.3e\n', var(w_after));
fprintf('ESS after update: %.2f  (N = %d)\n', pf.computeESS(), pf.N);

%% --- PLOT: the structure check ---
figure;
scatter(pf.particles(:,1), pf.particles(:,2), 25, w_after, 'filled');
colorbar; colormap(parula); hold on;
plot(x_s_true, y_s_true, 'rp', 'MarkerSize', 18, 'MarkerFaceColor', 'r');
plot(sensor.x,  sensor.y,  'ks', 'MarkerSize', 12, 'MarkerFaceColor', 'k');
xlabel('x'); ylabel('y');
title('Particle weights after one update (red star = true source, black square = sensor)');
axis([0 10 0 8]); axis equal; hold off;
%% --- RESAMPLE CHECK ---
fprintf('--- resample check ---\n');
ess_before = pf.computeESS();
nUnique_before = size(unique(pf.particles, 'rows'), 1);

pf.resample();

ess_after  = pf.computeESS();
nUnique_after = size(unique(pf.particles, 'rows'), 1);

fprintf('ESS     before / after : %.2f  / %.2f\n', ess_before, ess_after);
fprintf('unique  before / after : %d / %d\n', nUnique_before, nUnique_after);

%% --- ENTROPY CHECK ---
fprintf('--- entropy check ---\n');
fprintf('H prior  : %.4f   (log N = %.4f)\n', H_prior, log(pf.N));
fprintf('H post   : %.4f\n', H_post);
fprintf('delta H  : %.4f\n', H_prior - H_post);