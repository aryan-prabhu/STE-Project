function L2_masked = run_validation(dx, dy,boundary_type)
    % Runs the solver to steady state at a given grid spacing,
    % then returns the masked relative L2 error against the
    % analytical Gaussian plume (source neighbourhood excluded).

    % --- fixed physical domain; finer dx => more grid cells ---
    Lx = 40;  Ly = 60;
    Nx = round(Lx/dx);  Ny = round(Ly/dy);

    dt = 0.01;  u = 0.5;  v = 0.0;  D = 0.2;  Q = 10;
    src_x = 15;  src_y = 30;

    solver = AdvectionDiffusion2D(Nx,Ny,dx,dy,dt,u,v,D,src_x,src_y);
    solver.boundary_type = boundary_type;   % whatever this run wants
    % --- step forward until the field stops changing ---
    tol = 1e-6;
    for k = 1:200000
        C_old  = solver.C;
        solver = solver.step(Q);
        residual = max(abs(solver.C(:)-C_old(:))) / max(abs(solver.C(:)));
        if residual < tol
            break
        end
    end

    % --- analytical reference on the same grid ---
    src_i = solver.src_i;  src_j = solver.src_j;
    C_analytic = gaussianPlume(Nx,Ny,dx,dy,u,D,Q,src_i,src_j);
    C_num = solver.C;

    % --- mask out a fixed-physical-size neighbourhood of the source ---
    exclude_radius = 5/dx;
    [I,J] = ndgrid(1:Nx,1:Ny);
    dist_to_src = sqrt((I-src_i).^2 + (J-src_j).^2);
    mask = dist_to_src > exclude_radius;

    % --- masked relative L2 error (named vars, then indexed) ---
    err_field = (C_num - C_analytic) .* mask;
    ref_field =  C_analytic          .* mask;
    L2_masked = norm(err_field(:)) / norm(ref_field(:));
end