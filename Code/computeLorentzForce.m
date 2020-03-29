function [Fx, Fy] = computeLorentzForce(x_previous, y_previous)

global C;
global Efield_x;
global Efield_y;
global meshSize;
global nx;
global ny;

% -------- Maps electron positions to mesh positions -------- %
nx_electrons = ceil(x_previous/(meshSize*1e-9));
ny_electrons = ceil(y_previous/meshSize*1e-9);

Ex_electrons = zeros(1, length(x_previous));
Ey_electrons = zeros(1, length(y_previous));

% -------- Assign electric field vector to each electron -------- %
for i = 1:length(nx_electrons)
    
    if (nx_electrons(i) == 0) % this occationally happens for electrons whose positions have not yet been corrected by region bounndary conditions
        nx_electrons(i) = 1;
        
    elseif (nx_electrons(i) > nx)
        nx_electrons(i) = nx;
        
    end
    
    if (ny_electrons(i) == 0)
        ny_electrons(i) = 1;
        
    elseif (ny_electrons(i) > ny)
        ny_electrons(i) = ny;
        
    elseif (ny_electrons(i) < 0)
        ny_electrons(i) = 1;
        
    end
    
    Ex_electrons(i) = Efield_x(ny_electrons(i), nx_electrons(i));
    Ey_electrons(i) = Efield_y(ny_electrons(i), nx_electrons(i));
        
end

% -------- Compute force on all electrons and return -------- %
Fx = (C.q)*Ex_electrons;
Fy = (C.q)*Ey_electrons;

