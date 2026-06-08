% test_levelset.m — first check of LevelSetExtractor on a steady plume.
% --- build a small steady no-wind plume ---
Nx = 50; Ny = 50;
dx = 1;  dy = 1;
dt = 0.1;
u  = 6;  v = 3;          % no wind
D  = 0.5;                % isotropic diffusion
src_i = 25; src_j = 25;  % source near the centre

solver = AdvectionDiffusion2D(Nx, Ny, dx, dy, dt, u, v, D, src_i, src_j);

% --- run it forward to a near-steady field ---
Q = 1.0;
for k = 1:2000
    solver = solver.step(Q);
end

% --- extract a level set ---
C_th = 0.2 * max(solver.C(:));   % threshold at 20% of the peak value
extractor = LevelSetExtractor(dx, dy);
contours  = extractor.extract(solver.C, C_th);

% --- report ---
fprintf('Number of curves found: %d\n', numel(contours));
fprintf('Points in first curve:  %d\n', numel(contours(1).x));
fprintf('First point: (%.2f, %.2f)\n', contours(1).x(1),   contours(1).y(1));
fprintf('Last point:  (%.2f, %.2f)\n', contours(1).x(end), contours(1).y(end));

% --- plot field + contour ---
figure;
imagesc(solver.C'); axis equal tight; colorbar; hold on;
set(gca, 'YDir', 'normal');
for c = 1:numel(contours)
    plot(contours(c).x / dx + 1, contours(c).y / dy + 1, 'r-', 'LineWidth', 2);
end
title('Concentration field with extracted level set');
% --- quantitative check: is the contour really a circle? ---
src_x = (src_i-1) * dx;      % source physical x   (see note below)
src_y = (src_j-1) * dy;      % source physical y
radii = sqrt((contours(1).x - src_x).^2 + (contours(1).y - src_y).^2);
fprintf('Contour radius: mean %.3f, std %.4f, min %.3f, max %.3f\n', ...
        mean(radii), std(radii), min(radii), max(radii));