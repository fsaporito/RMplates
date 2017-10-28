function [gx,gy]=gradhphih2(i,x,y)
%
% gradienti delle funzioni di base
% sull'elemento di riferimento
% per k = 2
%
switch i
    case 1
        gx = 4*x + 4*y - 3;
        gy = 4*x + 4*y - 3;
    case 2
        gx = 4*x - 1;
        gy = 0;
    case 3
        gx = 0;
        gy = 4*y - 1;
    case 4
        gx = 4*y;
        gy = 4*x;
    case 5
        gx = -4*y;
        gy =  4 - 8*y - 4*x;
    case 6
        gx = 4 - 4*y - 8*x;
        gy = -4*x;
end
