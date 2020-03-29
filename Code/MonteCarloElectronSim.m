function [avgVelocityX, t_vec] = MonteCarloElectronSim(BottleNeck, DiffuseCollisions, ParticleScattering, FieldType, bottleNeckStartX, bottleNeckStartY, bottleneckWidth, bottleneckLength, doPlot)

% -------- Declare globals and define constants/initial parameters -------- %

global C;
global numElectrons;
global x_pos_init;
global y_pos_init;
global V_x_init;
global V_y_init;
global vx
global vy;
global x;
global y;
global boxWidthScaleFactor;
global boxLengthScaleFactor;
global maxTimeStep;
global Vapplied;

C.m_o = 9.10956e-31; %kg%
C.m = 0.26*C.m_o; %kg%
C.T = 300; %K%
C.k_b = 1.38064852e-23; %m^2 kg s^-2 K^-1%
C.q = 1.602e-19; %c%

meanTime = 0.2e-12; %s%
numElectrons = 1000;
maxTimeStep = 1000;
d = 100*boxWidthScaleFactor;

boxW = bottleneckWidth*100*boxWidthScaleFactor;
boxL = bottleneckLength*100*boxLengthScaleFactor;
verticalSep = (100*boxLengthScaleFactor) - (2*boxL);

% -------- Initialize positions and velocities of electrons -------- %

if (BottleNeck)
    initPositionAndVelocity("MB", "BottleNeck", bottleNeckStartX, bottleNeckStartY, bottleneckWidth, bottleneckLength);
else
    initPositionAndVelocity("MB", "Uniform", bottleNeckStartX, bottleNeckStartY, bottleneckWidth, bottleneckLength);
end

% -------- Set up more variables for computing electron trajectotries -------- %

avgVelocity_init = sum(sqrt((V_x_init.^2)+(V_y_init.^2)))/numElectrons; %average of maxwell boltzmann distribution%
dt = (1/500)*100*boxWidthScaleFactor/avgVelocity_init;

dx = zeros(1, numElectrons);
dy = zeros(1, numElectrons);
dvx = zeros(1, numElectrons);
dvy = zeros(1, numElectrons);
vx = zeros(1, numElectrons);
vy = zeros(1, numElectrons);
x = zeros(1, numElectrons);
y = zeros(1, numElectrons);

vx = V_x_init;
vy = V_y_init;
x = x_pos_init;
y = y_pos_init;
T = zeros(1, maxTimeStep);
t_vec = linspace(0, maxTimeStep*dt, maxTimeStep);

% randomly choose 10 electrons to have their trajectories over the whole
% simulation tracked
numElectronsToTrack = 10;
electronsToTrack = randi(numElectrons, 1, numElectronsToTrack);

x_oldvals = zeros(maxTimeStep, numElectronsToTrack);
y_oldvals = zeros(maxTimeStep, numElectronsToTrack);

%Calculate scattering probability%
Pscat = 1 - exp(-dt/meanTime); 

avgVelocity = zeros(1, maxTimeStep);
avgVelocityX = zeros(1, maxTimeStep);

% -------- Computation of particle acceleration for uniform electric field -------- %

if (FieldType == "Uniform")
    
    Ex = Vapplied/d; %V/m%
    Fx = (C.q)*Ex; %Magnitude of electric force - N% 
    ax = ones(1, numElectrons)*(Fx/(C.m)); %accceleration felt by electrons - m/s^2%
    ay = 0;
    
end

figure;

% -------- MAIN PROGRAM LOOP -------- %

for t = 1:maxTimeStep
    
% -------- Capturing previous electron positions -------- %

    x_previous = x;
    y_previous = y;

% -------- Update particle positions by solving Newton's equations (first order integration) -------- %
    
    % compute accelerations from finite difference E-field
    if (FieldType == "Custom")
        
        [Fx, Fy] = computeLorentzForce(x_previous, y_previous);
        ax = (Fx./(C.m))./(1e-9);
        ay = (Fy./(C.m))./(1e-9);
        
    end
    
    dvx = ax*dt;
    dvy = ay*dt;
    
    vx = vx + dvx;
    vy = vy + dvy;

    dx = (vx*dt) + ((ax*(dt^2))/2);
    dy = (vy*dt) + ((ay*(dt^2))/2);
    
    x = x_previous + dx;
    y = y_previous + dy;
    
    x_oldvals(t,:) = x_previous(electronsToTrack);
    y_oldvals(t,:) = y_previous(electronsToTrack);

% -------- Apply boundary conditions to applicable electrons -------- %
    
    toReflectY = find((y > 100*boxLengthScaleFactor) | (y < 0));
 
    toShiftRight = find(x < 0);
    toShiftLeft = find(x > 100*boxWidthScaleFactor);
    
    if ~isempty(toReflectY)
        vy(toReflectY) = -vy(toReflectY);
    end
    
    if ~isempty(toShiftRight)
        x(toShiftRight) = 100*boxWidthScaleFactor; 
    end
    
    if ~isempty(toShiftLeft)
        x(toShiftLeft) = 0;
       
    end
    
    % reflection off of bottleneck region walls   
    if (BottleNeck)
  
        for i = 1:numElectrons
            
            if (((x(i) >= bottleNeckStartX) && (x(i) <= (bottleNeckStartX + boxW))) && ((y(i) >= (bottleNeckStartY + boxL + verticalSep))))
                
                if (x_previous(i) <= bottleNeckStartX) && (y_previous(i) >= (bottleNeckStartY + boxL + verticalSep))
                    
                    if (DiffuseCollisions)
                        [vx(i), vy(i)] = thermalize(1);
                    else
                        vx(i) = -vx(i);
                    end
                    
                    x(i) = bottleNeckStartX;  
                    
                elseif (x_previous(i) >= (bottleNeckStartX + boxW)) && (y_previous(i) >= (bottleNeckStartY + boxL + verticalSep))
                    
                    if (DiffuseCollisions)
                        [vx(i), vy(i)] = thermalize(1);
                    else
                        vx(i) = -vx(i);
                    end
                    
                    x(i) = (bottleNeckStartX + boxW);
                               
                elseif (x_previous(i) > bottleNeckStartX) && (x_previous(i) < (bottleNeckStartX + boxW)) && (y_previous(i) >= (bottleNeckStartY + boxL + verticalSep))
                    
                    if (DiffuseCollisions)
                        [vx(i), vy(i)] = thermalize(1);
                    else
                        vy(i) = -vy(i);
                    end
                    
                    y(i) = (bottleNeckStartY + boxL + verticalSep);
                
                elseif (x_previous(i) > bottleNeckStartX) && (x_previous(i) < (bottleNeckStartX + boxW)) && (y_previous(i) < (bottleNeckStartY + boxL + verticalSep))
                    
                    if (DiffuseCollisions)
                        [vx(i), vy(i)] = thermalize(1);
                    else
                        vy(i) = -vy(i);
                    end
                    
                    y(i) = (bottleNeckStartY + boxL + verticalSep);
                    
                end
            
            elseif (((x(i) >= bottleNeckStartX) && (x(i) <= (bottleNeckStartX + boxW))) && ((y(i) <= (bottleNeckStartY + boxL))))
                
                if (x_previous(i) <= bottleNeckStartX) && (y_previous(i) <= (bottleNeckStartY + boxL))
                    
                    if (DiffuseCollisions)
                        [vx(i), vy(i)] = thermalize(1);
                    else
                        vx(i) = -vx(i);
                    end
                    
                    x(i) = bottleNeckStartX;
                                       
                elseif (x_previous(i) >= (bottleNeckStartX + boxW)) && (y_previous(i) <= (bottleNeckStartY + boxL))
                    
                    if (DiffuseCollisions)
                        [vx(i), vy(i)] = thermalize(1);
                    else
                        vx(i) = -vx(i);
                    end
                    
                    x(i) = (bottleNeckStartX + boxW);
                                
                elseif (x_previous(i) > bottleNeckStartX) && (x_previous(i) < (bottleNeckStartX + boxW)) && (y_previous(i) <= (bottleNeckStartY + boxL))
                    
                    if (DiffuseCollisions)
                        [vx(i), vy(i)] = thermalize(1);
                    else
                        vy(i) = -vy(i);
                    end
                    
                    y(i) = (bottleNeckStartY + boxL);
                    
                elseif (x_previous(i) > bottleNeckStartX) && (x_previous(i) < (bottleNeckStartX + boxW)) && (y_previous(i) > (bottleNeckStartY + boxL))
                    
                    if (DiffuseCollisions)
                        [vx(i), vy(i)] = thermalize(1);
                    else
                        vy(i) = -vy(i);
                    end
                    
                    y(i) = (bottleNeckStartY + boxL);
                    
                end
                
            end
             
        end  
    end
    
% -------- Random electron scattering -------- %. 
   
    % If any electrons
    % scatter, re-thermalize them (sample new speeds from maxwell boltzmann distribution)   
    
    if (ParticleScattering)
        
        scatteringProbabilites = rand(1, numElectrons);
        electronsToScatter = find(scatteringProbabilites <= Pscat);
        [vx(electronsToScatter), vy(electronsToScatter)] = thermalize(length(electronsToScatter));    
        
    end

% -------- calculate average electron velocities -------- %
    
    avgVelocity(t) = sum(sqrt((vx.^2)+(vy.^2)))/numElectrons; %average of maxwell boltzmann distribution%
    avgVelocityX(t) = (sum(vx))/numElectrons;
    
% -------- live plot of electrons moving around -------- %
    
    plot(x, y, 'b.');
    if (BottleNeck)
        hold on;
        plot([bottleNeckStartX bottleNeckStartX], [(bottleNeckStartY + boxL + verticalSep) 100*boxLengthScaleFactor], 'k');
        hold on;
        plot([(bottleNeckStartX + boxW) (bottleNeckStartX + boxW)], [(bottleNeckStartY + boxL + verticalSep) 100*boxLengthScaleFactor], 'k');
        hold on;
        plot([bottleNeckStartX bottleNeckStartX], [bottleNeckStartY (bottleNeckStartY + boxL)], 'k');
        hold on;
        plot([(bottleNeckStartX + boxW) (bottleNeckStartX + boxW)], [bottleNeckStartY (bottleNeckStartY + boxL)], 'k');
        hold on;
        plot([bottleNeckStartX (bottleNeckStartX + boxW)], [(bottleNeckStartY + boxL + verticalSep) (bottleNeckStartY + boxL + verticalSep)], 'k');
        hold on;
        plot([bottleNeckStartX (bottleNeckStartX + boxW)], [(bottleNeckStartY + boxL) (bottleNeckStartY + boxL)], 'k');
        hold off;
    end
    axis([0 100*boxWidthScaleFactor 0 100*boxLengthScaleFactor]);
    pause(0.0001);
    
end

% ------- Trajectory of select electrons over entire simulation duration ------- %
if (doPlot)
    figure;
    title('Particle Trajectories for 15 Randomly Selected Electrons');
    for pltCnt = 1:numElectronsToTrack

       x_diffs = diff(x_oldvals(:,pltCnt));
       x_diffs(length(x_diffs) + 1) = 0;

       rejections = find(abs(x_diffs) >= 100e-9);
       x_oldvals(rejections,pltCnt) = NaN;

       plot(x_oldvals(:,pltCnt), y_oldvals(:,pltCnt));
       hold on;
    end
    if (BottleNeck)
        plot([bottleNeckStartX bottleNeckStartX], [(bottleNeckStartY + boxL + verticalSep) 100*boxLengthScaleFactor], 'k');
        hold on;
        plot([(bottleNeckStartX + boxW) (bottleNeckStartX + boxW)], [(bottleNeckStartY + boxL + verticalSep) 100*boxLengthScaleFactor], 'k');
        hold on;
        plot([bottleNeckStartX bottleNeckStartX], [bottleNeckStartY (bottleNeckStartY + boxL)], 'k');
        hold on;
        plot([(bottleNeckStartX + boxW) (bottleNeckStartX + boxW)], [bottleNeckStartY (bottleNeckStartY + boxL)], 'k');
        hold on;
        plot([bottleNeckStartX (bottleNeckStartX + boxW)], [(bottleNeckStartY + boxL + verticalSep) (bottleNeckStartY + boxL + verticalSep)], 'k');
        hold on;
        plot([bottleNeckStartX (bottleNeckStartX + boxW)], [(bottleNeckStartY + boxL) (bottleNeckStartY + boxL)], 'k');      
    end
    xlabel('Horizontal Location (m)');
    ylabel('Vertical Location (m)');
    axis([0 100*boxWidthScaleFactor 0 100*boxLengthScaleFactor]);
    hold off;
end

% -------- Electron density map -------- %
finalPos = [transpose(x) transpose(y)];

if (doPlot)
    figure;
    hist3(finalPos, [100 100]);
    title('Electron Density Map at End of Simulation (Created Using 3D Histogram of Final Positions)');
    xlabel('x-Direction Bins (Bin Size = 2 nm)');
    ylabel('y-Direction Bins (Bin Size = 1 nm)');
    zlabel('Number of Electrons/bin');
end

% -------- Electron temperature map -------- %
finalTemps = ((vx.^2) + (vy.^2)).*((C.m)./(pi.*(C.k_b)));
xVals = linspace(0, 100*boxWidthScaleFactor, 200);
yVals = linspace(0, 100*boxLengthScaleFactor, 100);

[XVals, YVals] = meshgrid(xVals, yVals);
tempMap = griddata(x, y, finalTemps, XVals, YVals);

if (doPlot)
    figure;
    surf(XVals, YVals, tempMap);
    title('Temperature Map of Electrons at End of Simulation');
    xlabel('x-Position (nm)');
    ylabel('y-Position (nm)');
    zlabel('Temperature of Electrons');
end

% ------- Temperature plot over time ------- %

% T = ((avgVelocity.^2).*(C.m))./(pi.*(C.k_b)); 
% 
% figure;
% plot(t_vec, T, 'r');
% title("Semiconductor Temperature Vs Time");
% xlabel("Time (s)");
% ylabel("Temperature (K)")
% grid on;

