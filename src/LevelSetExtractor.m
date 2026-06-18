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
            %
            % FIXED (was swapped): contourc(Z, ...) always treats COLUMNS
            % of Z as its "X" and ROWS of Z as its "Y" -- same convention
            % as surf/mesh/imagesc. This project stores C(i,j) with i =
            % x-index (row) and j = y-index (column), per Equation 2 of
            % Report I (c_{i+1,j} increments i for the x-derivative). That
            % is the OPPOSITE of what contourc assumes, so contourc's
            % "raw row 1" is actually our column index (our y), and its
            % "raw row 2" is actually our row index (our x). The previous
            % version of this function used raw(1,:) for .x and raw(2,:)
            % for .y directly -- exactly backwards. This was invisible in
            % the Section 4 validation (centred source, no wind, isotropic
            % diffusion) because that field is symmetric under swapping
            % i<->j; it only shows up on an asymmetric (wind-driven) field.

            % --- call MATLAB's Marching Squares ---
            raw = contourc(C, [C_th C_th]);

            % --- unpack the packed matrix ---
            contours = struct('x', {}, 'y', {});
            col = 1;
            while col < size(raw, 2)
                n_pts = raw(2, col);                      % header: how many points follow
                idx_col = raw(1, col+1 : col+n_pts);       % contourc's "X" = column index of C = our j (y-index)
                idx_row = raw(2, col+1 : col+n_pts);       % contourc's "Y" = row index of C    = our i (x-index)

                % --- convert index space -> physical space (swap applied) ---
                phys_x = (idx_row - 1) * obj.dx;   % row index (our i, x) scaled by dx
                phys_y = (idx_col - 1) * obj.dy;   % column index (our j, y) scaled by dy

                contours(end+1) = struct('x', phys_x, 'y', phys_y);
                col = col + n_pts + 1;                     % jump to the next header
            end
        end
    end
end