classdef LevelSetExtractor
    properties
        dx
        dy
    end
    methods
        function obj = LevelSetExtractor(dx, dy)
            obj.dx = dx;
            obj.dy = dy;
        end

        function contours = extract(obj, C, C_th)
            % Extract the level set of field C at threshold C_th.
            % Returns a struct array; each element is one contour curve
            % with fields .x and .y in PHYSICAL coordinates (metres).

            % --- call MATLAB's Marching Squares ---
            raw = contourc(C, [C_th C_th]);

            % --- unpack the packed matrix ---
            contours = struct('x', {}, 'y', {});
            col = 1;
            while col < size(raw, 2)
                n_pts = raw(2, col);                  % header: how many points follow
                idx_x = raw(1, col+1 : col+n_pts);    % index-space x of those points
                idx_y = raw(2, col+1 : col+n_pts);    % index-space y of those points

                % --- convert index space -> physical space ---
                phys_x = (idx_x - 1) * obj.dx;
                phys_y = (idx_y - 1) * obj.dy;

                contours(end+1) = struct('x', phys_x, 'y', phys_y);
                col = col + n_pts + 1;                % jump to the next header
            end
        end
    end
end