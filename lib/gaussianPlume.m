function C_analytic = gaussianPlume(Nx, Ny, dx, dy, u, D, Q, src_i, src_j)
    % Analytical 2D Gaussian plume (far-field approximation).
    % Steady-state solution of the advection-diffusion PDE for a
    % continuous point source, uniform wind u in +x, isotropic D.

    C_analytic = zeros(Nx, Ny);

    for i = 1:Nx
        for j = 1:Ny
            % position of this cell relative to the source, in physical units
            x = (i - src_i) * dx;     % downwind distance
            y = (j - src_j) * dy;     % cross-wind distance

            if x > 0
                sigma2 = 2*D*x / u;                       % cross-wind variance
                amp    = (Q/u) * sqrt(u / (4*pi*D*x));    % centreline strength
                C_analytic(i,j) = amp * exp(-y^2 / (2*sigma2));
            else
                C_analytic(i,j) = 0;   % formula undefined at/upwind of source
            end
        end
    end
end