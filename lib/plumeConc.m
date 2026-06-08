function c = plumeConc(x_s, y_s, px, py, u, D, Q)
    % Point evaluation of the analytical 2D Gaussian plume.
    % Returns the steady-state concentration at ONE sensor point,
    % for a continuous point source at ONE location.
    %
    %   x_s, y_s : source position      (physical coordinates)
    %   px, py   : sensor position      (physical coordinates)
    %   u        : wind speed in +x
    %   D        : isotropic diffusion coefficient
    %   Q        : source strength
    %
    % Same physics as gaussianPlume.m, evaluated once instead of over a grid.

    x = px - x_s;          % downwind displacement, source -> sensor
    y = py - y_s;          % cross-wind displacement, source -> sensor

    if x > 0
        sigma2 = 2 * D * x / u;                 % cross-wind variance
        amp    = (Q / u) * sqrt(u / (4*pi*D*x)); % centreline strength
        c      = amp * exp( -y^2 / (2*sigma2) );
    else
        c = 0;             % sensor upwind of (or level with) source
    end
end