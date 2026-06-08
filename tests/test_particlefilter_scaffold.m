% test_particlefilter_scaffold.m
% Day 7 validation: confirm the ParticleFilterSIR scaffold builds a
% cloud that lies entirely inside the domain, with a valid weight set.

% --- domain bounds (must match the constructor call) ---
x_min = 0;  x_max = 10;
y_min = 0;  y_max = 8;
N = 1000;

% --- build the filter ---
pf = ParticleFilterSIR(N, x_min, x_max, y_min, y_max);

% --- check 1: every particle inside the domain box ---
xs = pf.particles(:, 1);
ys = pf.particles(:, 2);

x_ok = all(xs >= x_min) && all(xs <= x_max);
y_ok = all(ys >= y_min) && all(ys <= y_max);

fprintf('x range: [%.4f, %.4f]  (bounds [%.1f, %.1f])\n', ...
        min(xs), max(xs), x_min, x_max);
fprintf('y range: [%.4f, %.4f]  (bounds [%.1f, %.1f])\n', ...
        min(ys), max(ys), y_min, y_max);

% --- check 2: weights form a valid distribution ---
weight_sum   = sum(pf.weights);
weights_unif = all(pf.weights == pf.weights(1));

% --- verdict ---
if x_ok && y_ok && weights_unif && abs(weight_sum - 1) < 1e-12
    fprintf('PASS: cloud inside domain, weights uniform and sum to 1.\n');
else
    fprintf('FAIL: check the output above.\n');
end