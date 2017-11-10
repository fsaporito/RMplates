% ----------------------------------------------------
function [] = meshplot(xv,yv,endpoints,meshname,i);
% -----------------------------------------------------------------
% given some mesh data (see meshgen for format), it plots the mesh.
% -----------------------------------------------------------------
figure(i);
grid on;
title (meshname);
hold on

P = [xv,yv];
E = endpoints;
for j=1:size(E,1)
    x=[P(E(j,1),1) , P(E(j,2),1)];
    y=[P(E(j,1),2) , P(E(j,2),2)];
    plot(x,y)  
end

saveas (i, [meshname, '.png']);