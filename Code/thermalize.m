function [V_x,V_y] = thermalize(numElectrons)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    global C;

    vxProb = rand(1, numElectrons);
    vyProb = rand(1, numElectrons);
    
    mean = 0;
    s = sqrt((2*C.k_b)*(C.T)/(C.m));
    MB_dist = makedist('Normal', 'mu', mean, 'sigma', s);
    
    V_x = icdf(MB_dist, vxProb);
    V_y = icdf(MB_dist, vyProb);

end

