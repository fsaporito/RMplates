function [A, B, C] = RM_fem (...
                             xv, yv, ...
                             vertices, ...
                             edges, ...
                             endpoints, ...
                             fdq, ...
                             mu, ...
                             lambda, ...
                             t)
%
% -------------------------------------------------------------------------
% FEM solution of the Reissner-Midling plate bending
% Returns the three matrices that will compose the stiffness matrix
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
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% Quadrature Formula
% -------------------------------------------------------------------------

% (xhq,yhq) = quadrature node on reference element T_hat
% whq = relative quadrature nodes weight
[xhq yhq whq] = quadratura(fdq); % Quadrature Points Computation
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
    BE = zeros(12,6);     % Local B Matrix
    CE = zeros(6,6);      % Local C Matrix
    
    % A Matrix Computation Loop
    for i=1:12
        for j=1:12
            for q=1:Nq
                 if j < 6.5  % First six basis functions (j)
                   v_phi_j = [phihq(j,q); 0];
                   x_phi_j_1 = gphihqx(j,q);
                   y_phi_j_1 = gphihqy(j,q);
                   x_phi_j_2 = 0;
                   y_phi_j_2 = 0;
                elseif j > 6.5 % Second six basis functions (i);
                   v_phi_j = [0; phihq(j-6,q)];
                   x_phi_j_1 = 0;
                   y_phi_j_1 = 0;
                   x_phi_j_2 = gphihqx(j-6,q);
                   y_phi_j_2 = gphihqy(j-6,q);
                end
                if i < 6.5 % First six basis functions (i)
                    v_phi_i = [phihq(i,q); 0];
                    x_phi_i_1 = gphihqx(i,q);
                    y_phi_i_1 = gphihqy(i,q);
                    x_phi_i_2 = 0;
                    y_phi_i_2 = 0;
                elseif i > 6.5 % Second six basis functions (i)
                    v_phi_i = [0; phihq(i-6,q)];
                    x_phi_i_1 = 0;
                    y_phi_i_1 = 0;
                    x_phi_i_2 = gphihqx(i-6,q);
                    y_phi_i_2 = gphihqy(i-6,q);
                end
                
                % grad phi_j (it returns a 2x2 matrix)
                grad_phi_j_rif = JFIT*[x_phi_j_1 y_phi_j_1; 
                                       x_phi_j_2 y_phi_j_2];
                % grad phi_i (it returns a 2x2 matrix)
                grad_phi_i_rif = JFIT*[x_phi_i_1 y_phi_i_1; 
                                       x_phi_i_2 y_phi_i_2];
                
                % grad phi_i_t (it returns a 2x2 matrix)
                grad_phi_i_rif_t = [x_phi_i_1 x_phi_i_2; 
                                    y_phi_i_1 y_phi_i_2]*JFI;
                
                % A1 = div_ref_phi_j * div_ref_phi_i
                % div_ref_phi = trace(grad_phi)
                A1 = trace(grad_phi_j_rif) * trace(grad_phi_i_rif);
                
                % A2 = grad_phi_j : (grad_phi_i + grad_phi_i_transpost)
                A2 = contr(grad_phi_j_rif, ...
                           grad_phi_i_rif + grad_phi_i_rif_t);
                       
                % A3
                A3 = dot(v_phi_j,v_phi_i);
               
                % Full local A term
                AE(i,j) = AE(i,j) + (lambda*A1 + ...
                                     mu*A2 + ...
                                     mu*t^(-2)*A3)*whq(q);            
            end % End For q
            AE(i,j) = 2*area*AE(i,j);
        end % End For j
    end % End For i
    
    % B Matrix Computation Loop
    for i=1:12
        for k=1:6
            for q=1:Nq
                tmpJ = JFIT*[gphihqx(k,q); gphihqy(k,q)];
                if i < 6.5         % First six basis functions (i)
                    tmpI = [phihq(i,q); 0];
                elseif i > 6.5   % Second six basis functions (i)
                    tmpI = [0; phihq(i-6,q)];
                end
                BE(i,k) = dot(tmpJ,tmpI)*whq(q) + BE(i,k);           
            end % End For q
            BE(i,k) = 2*area*BE(i,k); % mu*t^(-2) after, multiplies whole B
        end % End For k
    end % End For i
    
     % C Matrix Computation Loop
    for l=1:6
        for k=1:6
            for q=1:Nq
                tmpK = JFIT*[gphihqx(k,q); gphihqy(k,q)];
                tmpL = JFIT*[gphihqx(l,q); gphihqy(l,q)];
                CE(l,k) = dot(tmpK,tmpL)*whq(q) + CE(l,k);           
            end % End For q
            CE(l,k) = 2*area*CE(l,k); % mu*t^(-2) after, multiplies whole C
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
    
    %disp(['dofg = ', num2str(dofg)]);
    %disp(['dofgg = ', num2str(dofgg)]);
    %disp(['sizeA = ', num2str(size(A))]);
    %disp(['sizeB = ', num2str(size(B))]);
    %disp(['sizeC = ', num2str(size(C))]);
    
    % Matrix A Asseblying
    A(dofgg,dofgg) = A(dofgg,dofgg) + AE;
    
    % Matrix B Asseblying
    B(dofgg,dofg) = B(dofgg,dofg) + BE;
    
    % Matrix C Asseblying
    C(dofg,dofg) = C(dofg,dofg) + CE;
    
end % End For

B = mu*t^(-2)*B;
C = mu*t^(-2)*C;


end % end function