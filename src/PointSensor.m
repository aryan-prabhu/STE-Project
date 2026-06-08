classdef PointSensor
    properties
        x       % sensor location, x-coordinate (metres)
        y       % sensor location, y-coordinate (metres)
        sigma   % noise standard deviation (instrument spec)
    end

    methods
        function obj = PointSensor(x, y, sigma)
            % constructor: store location and noise spec
            obj.x = x;
            obj.y = y;
            obj.sigma = sigma;
        end

        function c_true = sampleField(obj, solver)
            % bilinear interpolation of solver.C at (obj.x, obj.y)
            % returns the CLEAN concentration — no noise
            % --- locate the cell (0-based grid math) ---
            i = floor(obj.x / solver.dx);
            j = floor(obj.y / solver.dy);

            % --- fractional position inside the cell ---
            tx = (obj.x - i * solver.dx) / solver.dx;
            ty = (obj.y - j * solver.dy) / solver.dy;

            % --- convert to 1-based array indices ---
            ai = i + 1;
            aj = j + 1;

            % --- four corner values ---
            c00 = solver.C(ai,   aj);
            c10 = solver.C(ai+1, aj);
            c01 = solver.C(ai,   aj+1);
            c11 = solver.C(ai+1, aj+1);

            % --- bilinear blend ---
            c_bottom = (1 - tx) * c00 + tx * c10;
            c_top    = (1 - tx) * c01 + tx * c11;
            c_true   = (1 - ty) * c_bottom + ty * c_top;
        
        end

        function y_meas = measure(obj, solver)
            % calls sampleField, adds v_k ~ N(0, sigma^2)
            % returns the NOISY measurement
            c_true = obj.sampleField(solver);
            y_meas = c_true + obj.sigma * randn();
        end
    end
end