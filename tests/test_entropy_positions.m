clear classes

%% --- Truth forward model ---
Q   = 1.0;
u   = 1.0; v = 0.0; D = 0.5;
ad  = AdvectionDiffusion2D(50, 40, 1.0, 1.0, 0.1, u, v, D, 15, 20);
for k = 1:200
    ad = ad.step_loop(Q);
end

x = (0:ad.Nx-1) * ad.dx;
y = (0:ad.Ny-1) * ad.dy;

%% --- Plot truth field ---
figure;
imagesc(x, y, ad.C');
colorbar; axis xy;
title('Concentration field — truth');
xlabel('x (m)'); ylabel('y (m)');
hold on;
plot(15, 20, 'r*', 'MarkerSize', 12, 'LineWidth', 2);  % source
plot(20, 20, 'r*'); % inside
plot(30, 20, 'gs', 'MarkerSize', 12, 'LineWidth', 2);  % boundary
plot(45, 20, 'b^', 'MarkerSize', 12, 'LineWidth', 2);  % outside
legend('Source','Inside','Boundary','Outside');

%% --- Three sensor positions ---
sigma_noise = 0.01;
pos_A = [20, 20];  % on source
pos_B = [30, 20];  % mid-plume, centreline
pos_C = [45, 20];  % far downstream, centreline

positions  = {pos_A, pos_B, pos_C};
labels = {'Close to source(20,20)', 'Mid (30,20)', 'Far (45,20)'};

%% --- Run experiment ---
N = 1000;
H_prior_ref = log(N);
results = zeros(1, 3);

rng(42);  % seed for reproducibility

for trial = 1:3
    pos = positions{trial};

    % Fresh particle filter for each trial — same prior each time
    pf = ParticleFilterSIR(N, 0, ad.Nx*ad.dx, 0, ad.Ny*ad.dy);
    H_prior = pf.computeEntropy();

    % Sensor at this position
    sensor = PointSensor(pos(1), pos(2), sigma_noise);

    % One noisy measurement from the truth field
    y_k = sensor.measure(ad);

    % Weight update
    pf.update(y_k, sensor, u, D, Q);

    % Entropy BEFORE resampling
    H_post = pf.computeEntropy();
    results(trial) = H_prior - H_post;

    fprintf('Position %s: H_prior=%.4f, H_post=%.4f, DeltaH=%.4f\n', ...
        labels{trial}, H_prior, H_post, results(trial));
end

%% --- Bar chart ---
figure;
bar(results);
set(gca, 'XTickLabel', labels);
ylabel('Entropy reduction \DeltaH (nats)');
title('Information gain by UAV position — single measurement');
grid on;