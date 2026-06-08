classdef AdvectionDiffusion2D

    properties
        Nx
        Ny
        dx
        dy
        C
        t
        u
        v
        D
        dt
        src_i
        src_j
        boundary_type   % 'dirichlet' (zero walls) or 'absorbing' (zero-gradient)
        
    end

    methods
        function obj = AdvectionDiffusion2D(Nx, Ny, dx, dy, dt, u, v, D, src_i, src_j)
            % Constructor: assigns grid, physics, and source parameters.
            obj.Nx = Nx;
            obj.Ny = Ny;
            obj.dx = dx;
            obj.dy = dy;
            obj.dt = dt;
            obj.u  = u;
            obj.v  = v;
            obj.D  = D;

            % Snap source location to nearest grid index.
            obj.src_i = round(src_i/ dx);
            obj.src_j= round(src_j/ dy);

            % Initialise state: empty field, time zero.
            obj.C = zeros(Nx, Ny);
            obj.t = 0;
            obj.boundary_type = 'dirichlet';
            obj.checkStability();   % run the CFL check at construction
        end
        
        function obj = step_loop(obj, Q_value) %here for reference DO NOT DELETE not vectorised 
            % Advance the concentration field by one time step.
            % Q_value is the source strength to inject this step.

            C_old = obj.C;
            C_new = C_old;   % start as a copy; we overwrite the interior

            % Loop over interior cells only (skip the boundary row/column)
            for i = 2 : obj.Nx - 1
                for j = 2 : obj.Ny - 1

                    % --- Diffusion: central differences (2nd order) ---
                    diffusion = obj.D * ( ...
                        (C_old(i+1, j) - 2*C_old(i, j) + C_old(i-1, j)) / obj.dx^2 ...
                      + (C_old(i, j+1) - 2*C_old(i, j) + C_old(i, j-1)) / obj.dy^2 ...
                    );

                    % --- Advection: first-order upwind ---
                    if obj.u >= 0
                        adv_x = obj.u * (C_old(i, j) - C_old(i-1, j)) / obj.dx;
                    else
                        adv_x = obj.u * (C_old(i+1, j) - C_old(i, j)) / obj.dx;
                    end
                    if obj.v >= 0
                        adv_y = obj.v * (C_old(i, j) - C_old(i, j-1)) / obj.dy;
                    else
                        adv_y = obj.v * (C_old(i, j+1) - C_old(i, j)) / obj.dy;
                    end
                    advection = -(adv_x + adv_y);

                    % --- Forward Euler time update ---
                    C_new(i, j) = C_old(i, j) + obj.dt * (diffusion + advection);
                end
            end

            % --- Inject the source at its grid cell ---
            C_new(obj.src_i, obj.src_j) = C_new(obj.src_i, obj.src_j) ...
                                          + Q_value * obj.dt / (obj.dx * obj.dy);

            % Commit the new field and advance the clock
            obj.C = C_new;
            obj.t = obj.t + obj.dt;
        end
        
        function obj = step(obj, Q_value)
            C_old = obj.C;

            % ---- shifted sub-matrices for the stencil ----
            C_centre = C_old(2:end-1, 2:end-1);
            C_right  = C_old(3:end,   2:end-1);
            C_left   = C_old(1:end-2, 2:end-1);
            C_up     = C_old(2:end-1, 3:end);
            C_down   = C_old(2:end-1, 1:end-2);

            % ---- diffusion: central difference  ----
            diff_x = (C_right - 2*C_centre + C_left) / obj.dx^2;
            diff_y = (C_up    - 2*C_centre + C_down) / obj.dy^2;
            diffusion = obj.D * (diff_x + diff_y);

            % ---- advection: first-order upwind  ----
            if obj.u >= 0
                adv_x = obj.u * (C_centre - C_left)  / obj.dx;
            else
                adv_x = obj.u * (C_right  - C_centre) / obj.dx;
            end
            if obj.v >= 0
                adv_y = obj.v * (C_centre - C_down) / obj.dy;
            else
                adv_y = obj.v * (C_up     - C_centre) / obj.dy;
            end

            % ---- assemble new interior field ----
            C_new = C_old;
            C_new(2:end-1, 2:end-1) = C_centre ...
                + obj.dt * (diffusion - adv_x - adv_y);

            % ---- source injection (Day 2, unchanged) ----
            C_new(obj.src_i, obj.src_j) = C_new(obj.src_i, obj.src_j) ...
                + Q_value * obj.dt / (obj.dx * obj.dy);

            % ---- BOUNDARY CONDITION: applied after the stencil ----
            if strcmp(obj.boundary_type, 'absorbing')
                C_new(end, :) = C_new(end-1, :);   % right wall
                C_new(1,   :) = C_new(2,     :);   % left wall
                C_new(:, end) = C_new(:, end-1);   % top wall
                C_new(:, 1)   = C_new(:, 2);       % bottom wall
            elseif strcmp(obj.boundary_type, 'dirichlet')
                % edges stay at zero — nothing to do
            end

            obj.C = C_new;
            obj.t = obj.t + obj.dt;
        end
        
        
        function checkStability(obj)
            courant     = abs(obj.u)*obj.dt/obj.dx + abs(obj.v)*obj.dt/obj.dy;
            diff_number = obj.D*obj.dt*(1/obj.dx^2 + 1/obj.dy^2);

            if courant > 1
                warning('CFL: Courant number %.3f exceeds 1 — advection unstable.', courant);
            end
            if diff_number > 0.5
                warning('CFL: diffusion number %.3f exceeds 0.5 — diffusion unstable.', diff_number);
            end
        end
    end

end
