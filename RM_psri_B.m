function [A, B, C] = RM_psri_B (...
                             xv, yv, ...
                             vertices, ...
                             edges, ...
                             endpoints, ...
                             fdq, ...
                             mu, ...
                             lambda, ...
                             t, ...
                             alpha)
%
% ------------------------------------------------------------------------------
% FEM solution of the Reissner-Midling plate bending (PSRI)
% Returns the three matrices from the PSRI correction
% Uses [P2,P2]-elements
% xv: Array of x-coordinates of the mesh vertex points
% yv: Array of y-coordinates of the mesh vertex points
% vertices: Array of the mesh vertex points
% edges: Array of the mesh edges of every triangular element
% endpoints: Array of all pair of points that create an edge
% fdq: Quadrature order. es: 'degree=5'
% mu: plate material parameter, double
% lambda: plate material parameter, double
% t: plate thickness, double
% alpha: PSRI coefficient parameter
% ------------------------------------------------------------------------------


% ------------------------------------------------------------------------------
% Quadrature Formula
% ------------------------------------------------------------------------------

% (xhq,yhq) = quadrature node on reference element T_hat
% whq = relative quadrature nodes weight
[xhq, yhq, whq] = quadratura(fdq); % Quadrature Points Computation
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


% -------------------------------------------------------------------------
% Coefficient Matrices
% -------------------------------------------------------------------------

% Variable Definition
N = 2*(nver+nedge);
M = nver+nedge;
A = sparse(N,N);
B = sparse(N,M);
C = sparse(M,M);


% Computation For Loop
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
    
    % Triangle Barycenter
    xb = (x1+x2+x3)/3;
    yb = (y1+y2+y3)/3; 
    
    % Trasformation Map Jacobian
    JF = [x2-x1 x3-x1
          y2-y1 y3-y1];
          
    % Trasformation Map Jacobian Inverse
    JFI = inv(JF);
    
    % % Trasformation Map Jacobian Inverse Transposed
    JFIT = JFI';
    
    % Triangle Area
    area = 0.5*det(JF);
    
    % Local Variables Definition
    AE = zeros(12,12);    % Local A Matrix
    BE = zeros(12,6);    % Local B Matrix
    CE = zeros(6,6);    % Local C Matrix
    
    % A Matrix Computation Loop
    for i=1:12
        for j=1:12
            if j < 6.5  % First six basis functions (j)
                tmpJ = phih2(j,xb,yb);
             elseif j > 6.5 % Second six basis functions (i)
                tmpJ = phih2(j-6,xb,yb);
            end
            if i < 6.5  % First six basis functions (j)
                tmpI = phih2(i,xb,yb);
            elseif i > 6.5 % Second six basis functions (i)
                tmpI = phih2(i-6,xb,yb);
            end
            
            AE(i,j) = area*tmpJ*tmpI;         
            
        end % End For j
    end % End For i
    
    % B Matrix Computation Loop
    for i=1:12
        for k=1:6
            if i < 6.5  % First six basis functions (j)
                tmpI = [p_phi_1([v1 v2 v3],i); ...
                        0; 0; 0];
            elseif i > 6.5 % Second six basis functions (i)
                tmpI = [0; 0; 0; ...
                        p_phi_2([v1 v2 v3],i-6)];
            end
            BE(i,k) = (area/3)*dot(tmpI, [g_phi_1([v1 v2 v3],k); ...
                                          g_phi_2([v1 v2 v3],k)]);  
        end % End For
    end % End For i
    
     % C Matrix Computation Loop
    for l=1:6
        for k=1:6
           CE(l,k) = (area/3)*dot([g_phi_1([v1 v2 v3],l); ...
                                   g_phi_2([v1 v2 v3],l)], ...
                                   [g_phi_1([v1 v2 v3],k); ...
                                    g_phi_2([v1 v2 v3],k)]);
        end % End For k
    end % End For l
    
        
    % GLobal Assembling
    % Global degrees of freedom (dof):
    % - vertex i  ----> i
    % - edge l -----> nver+l
    
    % Triangle Edges
    l1 = edges(iele,1);
    l2 = edges(iele,2);
    l3 = edges(iele,3);
   
    % Global dof Array 
    dofg = [v1 v2 v3 nver+l1 nver+l2 nver+l3];   % one component
    dofgg = [dofg , dofg + (nver+nedge)];     % two components
    
    % Matrix A Asseblying
    A(dofgg,dofgg) = A(dofgg,dofgg) + AE;
    
    % Matrix B Asseblying
    B(dofgg,dofg) = B(dofgg,dofg) + BE;
    
    % Matrix C Asseblying
    C(dofg,dofg) = C(dofg,dofg) + CE;
    
end % End For

A = mu*(t^(-2) - alpha)*A;
B = mu*(t^(-2) - alpha)*B;
C = mu*(t^(-2) - alpha)*C;

end % end function
