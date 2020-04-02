function [V] = timeDomainNetworkSolver(timeStep, maxTimeSteps, inputSignal, Cmatrix, G, doPlot)
    
    global I_n;
    global question_5;

    V = zeros(7, maxTimeSteps);
    F = zeros(7, maxTimeSteps);
    
    F(1, :) = inputSignal;
    
    if question_5
        F(3, :) = -I_n;
        F(5, :) = -I_n;
    end
    
    if doPlot
        figure;
    end
    
    for t = 1:maxTimeSteps
        
        A = (Cmatrix/timeStep) + G;
        A_inv = inv(A);
        
        if t == 1
           V(:, t) = F(:, t); 
        
        else
        
            V(:, t) = A_inv*(Cmatrix*(V(:, t-1)/timeStep) + F(:, t));
            
        end
        
        
        if doPlot
            hold on;
            plot(linspace(timeStep, t*timeStep, (t*timeStep/timeStep)), V(1, 1:t), 'r');
            plot(linspace(timeStep, t*timeStep, (t*timeStep/timeStep)), V(5, 1:t), 'b');
            grid on;
            

            pause(0.0001);
        end
        
    end
    
    if doPlot
        legend('Input Signal', 'Output Voltage');
    end

end

