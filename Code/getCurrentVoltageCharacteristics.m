function [deviceResistance] = getCurrentVoltageCharacteristics()

global C;
global boxWidthScaleFactor;
global boxLengthScaleFactor;
global Efield_x;
global Efield_y;
global Vapplied;
global meshSize;

boxWidthScaleFactor = 2e-9; %nm%
boxLengthScaleFactor = 1e-9; %nm%

meshSize = 1; %nm%

chargesPerUnitArea = 1e19; % electrons/m^2;
Vsweep = [0.1 0.5 1 2 3 4 5 6 7 8 9 10];
avgVelocityX = zeros(1, length(Vsweep));

for i = 1:length(Vsweep)
    
    Vapplied = Vsweep(i);
    
    [Efield_x, Efield_y] = FiniteDifferenceSolver(1, 0.01, false, 80, 0, 1/5, 1.5/5);
    [avgVelocitiesX, dummy2] = MonteCarloElectronSim(true, false, true, "Custom", 40*boxWidthScaleFactor, 0, 1/5, 1.5/5, false);
    
    avgVelocityX(i) = sum(avgVelocitiesX)/length(avgVelocitiesX);
    
end

Jx = (C.q)*chargesPerUnitArea.*avgVelocityX;
Ix = Jx.*(100*boxLengthScaleFactor); % current in the x-direction over time

deviceFit = polyfit(Ix, Vsweep, 1);

figure;
plot(Ix, Vsweep, 'r.');
hold on;
plot(Ix, polyval(deviceFit, Ix), 'b');
grid on;
hold off;
xlim([0 max(Ix)]);
ylim([0 max(polyval(deviceFit, Ix))]);
title('Device Voltage vs. Current');
xlabel('Current (A)');
ylabel('Voltage (V)');
legend('Device Voltage vs. Current', 'Linear Fit to Data');

deviceResistance = deviceFit(1); % Ohms

end