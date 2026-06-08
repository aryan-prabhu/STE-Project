clear; clc;

C = [ 1  2  3  4 ;
      5  6  7  8 ;
      9 10 11 12 ;
     13 14 15 16 ];

[Nx, Ny] = size(C);

C_centre = C(2:Nx-1, 2:Ny-1);
C_right  = C(3:Nx,   2:Ny-1);
C_left   = C(1:Nx-2, 2:Ny-1);
C_up     = C(2:Nx-1, 3:Ny);
C_down   = C(2:Nx-1, 1:Ny-2);

x_curv = C_right - 2*C_centre + C_left
y_curv = C_up    - 2*C_centre + C_down

u  = -2.0;     % try +2.0 first, then rerun with -2.0
dx = 1.0;

if u >= 0
    dCdx = (C_centre - C_left ) / dx;
else
    dCdx = (C_right  - C_centre) / dx;
end

advection_x = -u * dCdx