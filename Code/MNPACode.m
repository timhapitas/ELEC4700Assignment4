global deviceResistance;

R1 = 1;
C = 0.25;
R2 = 2;
L = 0.2;
R3 = deviceResistance;
alpha = 100;
R4 = 0.1;
Ro = 1000;

Cmatrix = [0 0 0 0 0 0 0; C -C 0 0 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 -L 0];
G = [1 0 0 0 0 0 0; 1/R1 -((1/R1) + (1/R2)) 0 0 0 0 0; 0 0 1/R3 0 0 -1 0; 0 0 0 1 0 0 -alpha; 0 0 0 0 0 -1 1; 0 0 0 1/R4 -((1/R4) + (1/Ro)) 0 0; 0 1 -1 0 0 0 0];

% ------- DC Case ------- %

Vin = linspace(-10, 10, 20);
Vo = zeros(1, length(Vin));
V3 = zeros(1, length(Vin));

for i = 1:length(Vin)
   
    F = [Vin(i); 0; 0; 0; 0; 0; 0];
    
    V = G\F;
    
    Vo(i) = V(5);
    V3(i) = V(3);
    
end

figure;
plot(Vin, Vo, 'r');
hold on;
plot(Vin, V3, 'b');
hold off;
legend('Output Voltage', 'Node 3 Voltage');
grid on;
title('Plot of Output and Node 3 Voltages for a DC Sweep From -10 to 10 Volts');
xlabel('Input Voltage (V)');
ylabel('Solved Voltages (V)');

% ------- AC SWEEP ------- %

Vin = 10;
F = [Vin; 0; 0; 0; 0; 0; 0];
omega = linspace(0, 100, 1000);

Vo = zeros(1, length(omega));
V3 = zeros(1, length(Vin));

for count = 1:length(omega)
   
   tempMatrix = (G + 1j*omega(count)*Cmatrix);
   V = tempMatrix\F;
   
   Vo(count) = abs(V(5));
   V3(count) = abs(V(3));
    
end

figure;
plot(omega, Vo, 'r');
hold on;
plot(omega, V3, 'b');
hold off;
legend('Output Voltage', 'Node 3 Voltage');
grid on;
title('Plot of Output Voltage and Voltage at Node 3 Vs. Frequency');
xlabel('Frequency (rad/s)');
ylabel('Voltage (V)');

figure;
plot(omega, Vo/Vin, 'g');
grid on;
title('Plot of Voltage Gain Vs. Frequency');
xlabel('Frequency (rad/s)');
ylabel('Voltage Gain (V/V)');

% ------- Part d) ------- %

pdf = makedist('Normal', 'mu', C, 'sigma', 0.05);
omega = pi;
Vin = 10;
F = [Vin; 0; 0; 0; 0; 0; 0];

sampleCount = 1000;
C = zeros(1, sampleCount);
gain = zeros(1, sampleCount);

for i = 1:sampleCount
    
    C(i) = icdf(pdf, rand(1));
    Cmatrix = [0 0 0 0 0 0 0; -C(i) C(i) 0 0 0 0 0; 0 0 0 0 0 -L 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0];
    
    tempMatrix = (G + 1j*omega*Cmatrix);
    V = tempMatrix\F;
    
    gain(i) = abs(V(5))/Vin;

end

figure;
histogram(C, 15);
title('Histogram of Capacitance Values');
xlabel('Capacitance (F)');
ylabel('Number of Samples with Capacitance C');

figure;
histogram(gain, 15);
title('Histogram of Gain Values');
xlabel('Gain (V/V)');
ylabel('Number of Samples with Gain Vo/Vin');

% Re-assign variables for the remainder of the assignment
R1 = 1;
C = 0.25;
R2 = 2;
L = 0.2;
R3 = deviceResistance;
alpha = 100;
R4 = 0.1;
Ro = 1000;
