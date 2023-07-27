function [BRmap] = blueRedMap(n,grayVal)
%BRMAP makes a blue-gray-red colormap
%   BRMap(N) makes a 3 X N colormap ranging from blue to gray to red.
%       The grayVal argument defines the saturation of the gray in RGB
%       values from 0 to 255. 
%

% authordate: EHS 20160106

if mod(n,2)
    BRmap = [linspace(255,0,n)' [linspace(0,grayVal-1,(
    n-1)/2) grayVal linspace(grayVal-1,0,n-1/2)]' linspace(0,255,n)'];
else
    BRmap = [linspace(255,0,n)' [linspace(0,grayVal,n/2)  linspace(grayVal,0,n/2)]' linspace(0,255,n)'];
end

BRmap = BRmap./255;

end

