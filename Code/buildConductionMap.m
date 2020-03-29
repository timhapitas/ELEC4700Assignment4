function [cMap] = buildConductionMap(L, W, backgroundCond, boxCond, boxXstart, boxYstart, bottleneckLength, bottleneckWidth)

%"bottleneckLength" and "bottleneckWidth" are defined as fractions of the total
%region length and width.

global meshSize;
global nx;
global ny;

cMap = zeros(ny, nx);
boxL = bottleneckLength*L;
boxH = bottleneckWidth*W;
boxXfinish = boxXstart + boxL/meshSize;
boxYfinish = boxYstart + (boxH/meshSize);
verticalSep = W - (2*boxH/meshSize);

for j = 1:ny
    for i = 1:nx
       
        if (((i >= boxXstart) && (i <= boxXfinish)) && ((j >= boxYstart) && (j <= boxYfinish)))
            cMap(j,i) = boxCond;
            
        elseif (((i >= boxXstart) && (i <= boxXfinish)) && ((j >= (boxYfinish + verticalSep)) && (j <= ny)))
            cMap(j,i) = boxCond;
            
        else
            cMap(j,i) = backgroundCond;
            
        end    
    end
end
end

