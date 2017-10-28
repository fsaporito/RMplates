function [p_phi_1, p_phi_2, g_phi_1, g_phi_2] = L2proj_phi( ...
                                                            xv, yv, ...
                                                            vertices, ...
                                                            i)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nver  = length(xv);         % vertex number
nele = size(vertices,1);    % elements number
Ah = sparse(nver, nver); % Stiffness Matrix (later it will be duplicated)
bh = zeros(nver, 1); % Load Term phi (x and y components)
bhgx = zeros(nver, 1); % Load Term grad psi (x component)
bhgy = zeros(nver, 1); % Load Term grad psi (y component)

% Main Computation Loop On Every Triangle
for iele=1:nele
    
    % Vertex 1
    v1 = vertices(iele,1);
    x1 = xv(v1);
    y1 = yv(v1);
    
    % Vertex 2
    v2 = vertices(iele,2);
    x2 = xv(v2);
    y2 = yv(v2);
    
    % Vertex 3
    v3 = vertices(iele,3);
    x3 = xv(v3);
    y3 = yv(v3);
    
    % Triangle Area With The Determinant Formula
    T = 0.5*det([1  1  1
                 x1 x2 x3
                 y1 y2 y3]);
                      
    % KhT Element Matrix 
    KhT = (T/3)*[1 0 0; ...
                 0 1 0; ...
                 0 0 1];
                
    % Ah Matrix Computation
    Ah([v1 v2 v3], [v1 v2 v3]) = Ah([v1 v2 v3], [v1 v2 v3]) + KhT;
    
    % Elementary Load Term Computation With Trapezoid Method
    fhT_t = (T/3)*[phih2(i,x1,y1) phih2(i,x2,y2) phih2(i,x3,y3)]';
    [gx1 gy1] = gradhphih2(i,x1,y1);
    [gx2 gy2] = gradhphih2(i,x2,y2);
    [gx3 gy3] = gradhphih2(i,x3,y3);
    fhT_tgx = (T/3)*[gx1 gx2 gx3]';
    fhT_tgy = (T/3)*[gy1 gy2 gy3]';
    
    % General Constant Term fh Computation
    bh([v1 v2 v3]) = bh([v1 v2 v3]) + fhT_t;     
    bhgx([v1 v2 v3]) = bhgx([v1 v2 v3]) + fhT_tgx;
    bhgy([v1 v2 v3]) = bhgy([v1 v2 v3]) + fhT_tgy;
    
end

% Duplication (it's vectorial)
Kh = [Ah, zeros(nver,nver); ...
      zeros(nver,nver), Ah];
Fh = [bh; bh];
Fhg = [bhgx; bhgy];

% Solutions Array
uh = zeros(2*nver,1);
uhg = zeros(2*nver,1);

% Solution Computation
uh = Kh\Fh;
uhg = Kh\Fhg;

% Return The Coefficients Of The Two Components
p_phi_1 = uh(1:nver);
p_phi_2 = uh(nver+1,end);
g_phi_1 = uhg(1:nver);
g_phi_2 = uhg(nver+1:end);


end

