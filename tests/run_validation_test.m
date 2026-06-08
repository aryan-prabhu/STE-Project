clear; close all;

dx_values = [1, 0.5, 0.25];

for dx = dx_values
    e_dir = run_validation(dx, dx, 'dirichlet');
    e_abs = run_validation(dx, dx, 'absorbing');
    fprintf('dx = %.3f : dirichlet = %.2f%%  |  absorbing = %.2f%%\n', ...
            dx, 100*e_dir, 100*e_abs);
end