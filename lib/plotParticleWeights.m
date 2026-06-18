function plotParticleWeights(particles, weights, mode, opts)
%PLOTPARTICLEWEIGHTS Visualise a particle cloud weighted by posterior weight.
%
%   plotParticleWeights(particles, weights, mode)
%   plotParticleWeights(particles, weights, mode, opts)
%
%   particles : N x 2 array, columns [x_s, y_s]
%   weights   : N x 1 vector of particle weights (need not be pre-normalised)
%   mode      : 'color' | 'sqrt'
%
%       'color' - existing behaviour. Constant marker size; weight encoded
%                 by colour via the current colormap (linear scale, as in
%                 Figure 7 of Report I).
%
%       'sqrt'  - constant colour; marker AREA proportional to sqrt(weight).
%                 Area, not radius, is set proportional to sqrt(w) so that
%                 the *visual area* a viewer perceives scales linearly with
%                 weight (perceptual honesty for area-based encodings).
%
%   A 'log' mode was tried and removed after validation (single-step case,
%   Report I Section 5.5): log-compression treats every decade of weight as
%   an equally-sized step on the size axis, which gives the near-dead tail
%   (many decades, scientifically meaningless) as much size-range as the
%   one decade that actually distinguishes ridge from near-ridge. Result
%   was a visually blurred ridge boundary -- confirmed against 'sqrt',
%   which collapses the dead tail toward the minimum size correctly because
%   sqrt of a tiny number is still tiny in linear terms.
%
%   opts (optional struct), fields:
%       .truePos    [x_s, y_s] true source location -> plotted as a red star
%       .sensorPos  [x, y]     sensor location -> plotted as a black square
%       .sizeRange  [minPts, maxPts] marker size range in points^2
%                   (default [4, 200])
%       .title      string for the plot title
%
%   Validate against the single-step case of Report I, Section 5.5
%   (true source (3,4), sensor (6,4), N = 1000, uniform prior) BEFORE
%   reusing this in any multi-step experiment -- see
%   test_plotParticleWeights.m.

    if nargin < 4
        opts = struct();
    end
    if ~isfield(opts, 'sizeRange')
        opts.sizeRange = [4, 200];
    end

    w = weights(:);
    wSum = sum(w);
    if wSum <= 0
        error('plotParticleWeights:badWeights', ...
            'Weights sum to zero or are empty -- nothing to plot.');
    end
    w = w / wSum;   % defensive renormalisation; caller's weights should
                    % already sum to 1, but this guards against accidentally
                    % passing raw tilde-weights (Equation eq:wtilde) by mistake

    switch mode
        case 'color'
            sz   = 30 * ones(size(w));   % constant marker size
            cVal = w;                    % colour encodes weight

        case 'sqrt'
            sz   = mapToSizeRange(sqrt(w), opts.sizeRange);
            cVal = [];                    % single colour; size encodes weight

        otherwise
            error('plotParticleWeights:badMode', ...
                'mode must be ''color'' or ''sqrt''.');
    end

    figure; hold on;
    if strcmp(mode, 'color')
        scatter(particles(:,1), particles(:,2), sz, cVal, 'filled');
        colorbar;
        colormap(parula);
    else
        scatter(particles(:,1), particles(:,2), sz, [0.20 0.45 0.85], ...
            'filled', 'MarkerFaceAlpha', 0.65);
    end

    if isfield(opts, 'truePos')
        plot(opts.truePos(1), opts.truePos(2), 'r*', ...
            'MarkerSize', 14, 'LineWidth', 2);
    end
    if isfield(opts, 'sensorPos')
        plot(opts.sensorPos(1), opts.sensorPos(2), 'ks', ...
            'MarkerSize', 10, 'MarkerFaceColor', 'k');
    end

    axis equal; box on;
    if isfield(opts, 'title')
        title(opts.title);
    end
    hold off;
end

function sz = mapToSizeRange(v, sizeRange)
    vMin = min(v);
    vMax = max(v);
    if (vMax - vMin) < eps
        % degenerate case: all values equal (e.g. the uniform prior) --
        % every particle gets the same mid-range size rather than dividing
        % by zero.
        sz = mean(sizeRange) * ones(size(v));
    else
        frac = (v - vMin) / (vMax - vMin);
        sz = sizeRange(1) + frac * (sizeRange(2) - sizeRange(1));
    end
end