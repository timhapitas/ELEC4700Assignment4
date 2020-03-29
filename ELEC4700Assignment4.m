close all;
clear;
clc;

addpath('./code/');

global deviceResistance;

% -------- Question 1 and 2 -------- %
% deviceResistance = getCurrentVoltageCharacteristics();
deviceResistance = 40.6713; % Ohms - manually setting this equal to the result of getCurrentVoltageCharacterisitcs to avoid running the simulation
                            % assignment 3 over and over again.

% -------- Question 3 -------- %
MNPACode;

% -------- Question 4 -------- %

% circuit matrices
Cmatrix = [0 0 0 0 0 0 0; C -C 0 0 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 -L 0];
G = [1 0 0 0 0 0 0; 1/R1 -((1/R1) + (1/R2)) 0 0 0 0 0; 0 0 1/R3 0 0 -1 0; 0 0 0 1 0 0 -alpha; 0 0 0 0 0 -1 1; 0 0 0 1/R4 -((1/R4) + (1/Ro)) 0 0; 0 1 -1 0 0 0 0];

% some simulation parameters
simulationDuration = 1; %s%
maxTimeSteps = 1000;
timeStep = simulationDuration/maxTimeSteps;
samplingFrequency = 1/timeStep;
t_vec = 0:timeStep:(simulationDuration - timeStep);
freqRange = (0:maxTimeSteps-1)*(samplingFrequency/maxTimeSteps);

% defining input signals to feed into the time domain circuit solver
stepFunction = zeros(1, maxTimeSteps);
stepFunction(1, 30:end) = 1; %V%
harmonicInput = sin(2*pi*(1/0.03).*t_vec);
gaussianPulse = exp(-0.5*((t_vec - 0.06)./0.03).^2);

simSignals = {stepFunction, harmonicInput, gaussianPulse};
nameVec = ["Step Input Function", "Harmonic Input Function", "Gaussian Input Pulse"];

% plots of network response to the above defined input signals
for i = 1:max(size(simSignals))
    
    systemVariables = timeDomainNetworkSolver(timeStep, maxTimeSteps, simSignals{1, i}, Cmatrix, G, false);
    
    % perform fourier transform of time domain simulation signals
    inputSignalFT = abs(fftshift(fft(simSignals{1, i})));
    outputVoltageFT = abs(fftshift(fft(systemVariables(5,:))));
 
    figure;
    hold on;
    plot(t_vec, systemVariables(1,:), 'r');
    plot(t_vec, systemVariables(5,:), 'b');
    grid on;
    hold off;
    title('Time Domain Simulation Result - Output Voltage for a ' + nameVec(i));
    xlabel('Time Elapsed (s)');
    ylabel('Signal Voltage (V)');
    legend('Input Signal', 'Output Voltage');
    
    figure;
    hold on;
    plot((freqRange - maxTimeSteps/2), inputSignalFT, 'r');
    plot((freqRange - maxTimeSteps/2), outputVoltageFT, 'b');
    grid on;
    hold off;
    title('Fourier Components of Input and Output Signals');
    xlabel('Frequency (f)');
    ylabel('Magnitude of Fourier Transform');
    xlim([-100 100]);
    legend('Input Signal', 'Output Voltage');

end








