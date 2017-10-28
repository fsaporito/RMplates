function b = P2load(xv,yv,vertices,edges,boundary,boundedges,endpoints);
%
% Compute The Load Term

% ------------------------------------------------------------------------------
% Quadrature Formula
% ------------------------------------------------------------------------------

% (xhq,yhq) = quadrature node on reference element T_hat
% whq = relative quadrature nodes weight
fdq = 'degree=5'; % Quadrature Degree
[xhq,yhq,whq]=quadratura(fdq); % Quadrature Points Computation
Nq = length(xhq); % Number of quadrature points

% Mesh info
nver  = length(xv);         % vertex number
nedge = size(endpoints,1);  % edges number
nele = size(vertices,1);    % elements number

% Basis functions 
% computed at the quadrature points of the reference element
phihq = zeros(6,Nq);
for i=1:6 % Computation Loop
    for q=1:Nq
        phihq(i,q) = phih2(i,xhq(q),yhq(q));
    end % end for
end % end for

% Basis functions gradient 
% computed at the quadrature points of the reference element
gphihqx = zeros(6,Nq);
gphihqy = zeros(6,Nq);
for i=1:6 % Computation Loop
    for q=1:Nq
        [gx,gy] = gradhphih2(i,xhq(q),yhq(q));
        gphihqx(i,q) = gx;
        gphihqy(i,q) = gy;
    end % end for
end % end for

% ----------------

b = zeros(nver+nedge,1); 

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
    
    % Trasformation Map Jacobian
    JF = [x2-x1 x3-x1
          y2-y1 y3-y1];
          
    % Trasformation Map Jacobian Inverse
    JFI = inv(JF);
    
    % % Trasformation Map Jacobian Inverse Transposed
    JFIT = JFI';
    
    % Triangle Area
    area = 0.5*det(JF);
    
    % Local Load Term
    fE = zeros(6,1);
        
    for i=1:6
        for q=1:Nq
            %
            % Compute the image in (xq,yq) sul triangolo
            % corrente del nodo di quadratura
            % (xhq(q),yhq(q)) che sta
            % sull'elemento di rif.
            %
            tmp = JF*[xhq(q);yhq(q)]+[x1;y1];
            %
            xq = tmp(1);
            yq = tmp(2);
            %
            fE(i) = fE(i) + load_plate(xq,yq)*phihq(i,q)*whq(q);
        end % End For q
        %
        fE(i) = 2*area*fE(i);
        %
    end % End For i
    
    % GLobal Assembling
    % Global degrees of freedom (dof):
    % - vertex i  ----> i
    % - edge l -----> nver+l
    
    % Triangle Edges
    l1 = edges(iele,1);
    l2 = edges(iele,2);
    l3 = edges(iele,3);
   
    % Global dof Array     
    dofg = [v1 v2 v3 nver+l1 nver+l2 nver+l3];
    
    % Global Load Term Update
    b(dofg) = b(dofg) + fE;
   
end % End Function
