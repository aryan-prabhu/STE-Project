classdef ParticleFilterSIR < handle
    % Particle filter (SIR) for source localisation - Eq (4).
    % State is the particle cloud over theta = (x_s, y_s, Q) and the
    % weights. x_s, y_s drawn uniformly from the domain; Q drawn from
    % an informative Gamma(k=2, theta_Q) prior biased toward weak
    % releases (Hutchinson et al.-style scaled release-strength prior).

    properties
        N          % number of particles
        particles  % N-by-3 array: col 1 = x_s, col 2 = y_s, col 3 = Q
        weights    % N-by-1 array: one weight per particle
    end

    methods
        function obj = ParticleFilterSIR(N, x_min, x_max, y_min, y_max, theta_Q)
            if nargin < 6
                theta_Q = 0.25;
            end

            obj.N = N;

            xs = x_min + (x_max - x_min) * rand(N, 1);
            ys = y_min + (y_max - y_min) * rand(N, 1);

            U = rand(N, 2);
            Qs = -theta_Q * log(U(:,1) .* U(:,2));
            %Qs = gamrnd(2, theta_Q, N, 1);

            obj.particles = [xs, ys, Qs];
            obj.weights = ones(N, 1) / N;
        end

        function update(obj, y_k, sensor, u, D)
            g = zeros(obj.N, 1);

            for i = 1:obj.N
                x_s = obj.particles(i, 1);
                y_s = obj.particles(i, 2);
                Q_i = obj.particles(i, 3);

                c_pred = plumeConc(x_s, y_s, sensor.x, sensor.y, u, D, Q_i);

                r = y_k - c_pred;
                g(i) = exp( -r^2 / (2 * sensor.sigma^2) );
            end

            obj.weights = obj.weights .* g;
            obj.weights = obj.weights / sum(obj.weights);
        end

        function ess = computeESS(obj)
            ess = 1 / sum(obj.weights.^2);
        end

        function resample(obj)
            N = obj.N;
            edges = cumsum(obj.weights);
            edges(end) = 1;

            u0 = rand / N;
            pointers = u0 + (0:N-1)' / N;

            idx = zeros(N, 1);
            j = 1;
            for i = 1:N
                while pointers(i) > edges(j)
                    j = j + 1;
                end
                idx(i) = j;
            end

            obj.particles = obj.particles(idx, :);
            obj.weights   = ones(N, 1) / N;
        end

        function roughen(obj, K)
            if nargin < 2
                K = 0.2;
            end
            d = 3;
            Ex = max(obj.particles(:,1)) - min(obj.particles(:,1));
            Ey = max(obj.particles(:,2)) - min(obj.particles(:,2));
            EQ = max(obj.particles(:,3)) - min(obj.particles(:,3));
            sigma_x = K * Ex * obj.N^(-1/d);
            sigma_y = K * Ey * obj.N^(-1/d);
            sigma_Q = K * EQ * obj.N^(-1/d);
            obj.particles(:,1) = obj.particles(:,1) + sigma_x * randn(obj.N, 1);
            obj.particles(:,2) = obj.particles(:,2) + sigma_y * randn(obj.N, 1);
            obj.particles(:,3) = obj.particles(:,3) + sigma_Q * randn(obj.N, 1);
            obj.particles(:,3) = max(obj.particles(:,3), 0);
        end

        function H = computeEntropy(obj)
            w = obj.weights;
            w = w(w > 0);
            H = -sum(w .* log(w));
        end
    end
end