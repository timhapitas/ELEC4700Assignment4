function [] = initPositionAndVelocity(randnOrMB, UniformOrBottleNeck, bottleNeckStartX, bottleNeckStartY, bottleneckWidth, bottleneckLength)

global numElectrons;
global x_pos_init;
global y_pos_init;
global V_x_init;
global V_y_init;
global magThermalSpeed;
global boxWidthScaleFactor;
global boxLengthScaleFactor;

boxW = bottleneckWidth*100*boxWidthScaleFactor;
boxL = bottleneckLength*100*boxLengthScaleFactor;
verticalSep = (100*boxLengthScaleFactor) - (2*boxL);

if (UniformOrBottleNeck == "Uniform")
    
    x_pos_init = rand(1, numElectrons)*100*boxWidthScaleFactor;
    y_pos_init = rand(1, numElectrons)*100*boxLengthScaleFactor;

elseif (UniformOrBottleNeck == "BottleNeck")
    
    x_init = rand(1, numElectrons)*100*boxWidthScaleFactor;
    y_init = rand(1, numElectrons)*100*boxLengthScaleFactor;
    
    for i = 1:length(x_init)
    
        while ((x_init(i) >= bottleNeckStartX) && (x_init(i) <= (bottleNeckStartX + boxW)) && (y_init(i) >= (bottleNeckStartY + boxL + verticalSep)) || (x_init(i) >= bottleNeckStartX) && (x_init(i) <= (bottleNeckStartX + boxW)) && (y_init(i) <= (bottleNeckStartY + boxL)))

            x_init(i) = rand*100*boxWidthScaleFactor;
            y_init(i) = rand*100*boxLengthScaleFactor;

        end
        
    end
    
    x_pos_init = x_init;
    y_pos_init = y_init;

end

if(randnOrMB == "randn")

    directionAngles = randn(1, numElectrons)*2*pi;

    V_x_init = magThermalSpeed*cos(directionAngles);
    V_y_init = magThermalSpeed*sin(directionAngles);

elseif(randnOrMB == "MB")   
    
    [V_x_init, V_y_init] = thermalize(numElectrons);
    
else
    
    directionAngles = rand(1, numElectrons)*2*pi;

    V_x_init = magThermalSpeed*cos(directionAngles);
    V_y_init = magThermalSpeed*sin(directionAngles);
    
end

end