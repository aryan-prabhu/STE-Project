classdef ParticleFilterSIR < handle
    % Particle filter (SIR) for source localisation - Eq (4).
    % Day 7 scope: scaffold only. State is the particle cloud over
    % (x_s, y_s) and the weights. N particles drawn uniformly from
    % the prior; weights initialised uniform. Q held aside for now.

    properties
        N           % number of particles
        particles   % N-by-2 array: column 1 = x_s, column 2 = y_s
        weights     % N-by-1 array: one weight per particle
    end

    methods
        function obj = ParticleFilterSIR(N, x_min, x_max, y_min, y_max)
            obj.N = N;

            % Draw N particles uniformly at random from the domain.
            xs = x_min + (x_max - x_min) * rand(N, 1);
            ys = y_min + (y_max - y_min) * rand(N, 1);
            obj.particles = [xs, ys];

            % Uniform prior weights: every particle equally plausible.
            obj.weights = ones(N, 1) / N;
        end
        function update(obj, y_k, sensor, u, D, Q)
            % Eq (4): single measurement update — reweight particles by likelihood.
            % y_k    : the noisy sensor reading (scalar)
            % sensor : a PointSensor object (knows its own x, y, sigma)
            % u,v,D,Q: forward-model parameters, known constants for Day 8

            g = zeros(obj.N, 1);              % likelihood for each particle

            for i = 1:obj.N
                x_s = obj.particles(i, 1);    % this particle's hypothesised source x
                y_s = obj.particles(i, 2);    % this particle's hypothesised source y

                % predicted concentration at the sensor IF the source were here
                c_pred = plumeConc(x_s, y_s, sensor.x, sensor.y, u, D, Q);

                r = y_k - c_pred;                         % residual
                g(i) = exp( -r^2 / (2 * sensor.sigma^2) );% likelihood (taper)
            end

            obj.weights = obj.weights .* g;               % Bayes: prior x evidence
            obj.weights = obj.weights / sum(obj.weights); % renormalise to sum 1
        end
        function ess = computeESS(obj)
            % Effective Sample Size of the current weight vector.
            % ess in [1, N]: ~N means healthy, ~1 means fully degenerate.
            ess = 1 / sum(obj.weights.^2);
        end
        function resample(obj)
            % Systematic resampling. Draws N new particles with
            % replacement, prob. proportional to weight. Resets
            % weights to 1/N. Handle class: mutates obj in place.
            N = obj.N;

            % 1. Cumulative sum of weights: a [0,1] number line
            %    partitioned into N bins, bin i sized by weight i.
            edges = cumsum(obj.weights);
            edges(end) = 1;                 % guard against rounding

            % 2. One random start in [0, 1/N), then N evenly
            %    spaced pointers 1/N apart.
            u0 = rand / N;
            pointers = u0 + (0:N-1)' / N;

            % 3. For each pointer, find which bin it lands in.
            idx = zeros(N, 1);
            j = 1;
            for i = 1:N
                while pointers(i) > edges(j)
                    j = j + 1;
                end
                idx(i) = j;
            end

            % 4. Keep the chosen particles; reset weights.
            obj.particles = obj.particles(idx, :);
            obj.weights   = ones(N, 1) / N;
        end
        function H = computeEntropy(obj)
            % Shannon entropy of the weight vector.
            % H in [0, log(N)]: 0 = certain, log(N) = maximally uncertain.
            w = obj.weights;
            w = w(w > 0);                  % drop zeros: 0*log(0) := 0
            H = -sum(w .* log(w));
        end
    end
end