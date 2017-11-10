function [errL2, errH1] = RM (...
                               xv, yv, ...
                               vertices, ...
                               edges, ...
                               endpoints, ...
                               boundary, ...
                               boundedges, ...
                               fdq, ...
                               mu, ...
                               lambda, ...
                               t, ...
                               psri, ...
                               h, ...
                               out ...
                               )
%
% -------------------------------------------------------------------------
% FEM solution of the Reissner-Midling plate bending
% Uses [P2,P2]-elements
% xv: Array of x-coordinates of the mesh vertex points
% yv: Array of y-coordinates of the mesh vertex points
% vertices: Array of the mesh vertex points
% edges: Array of the mesh edges of every triangular element
% endpoints: Array of all pair of points that create an edge
% boundary: Points at the domain boundary
% boundedges: Edges at the domain boundary
% fdq: Quadrature order. es: 'degree=5'
% mu: plate material parameter, double
% lambda: plate material parameter, double
% t: plate thickness, double
% psri: flag, yes: use psri 
% out: flag, yes: verbose command line output, no: suppress output
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% Mesh info
% -------------------------------------------------------------------------

% Variable Definition
nver  = length(xv);         % vertex number
nedge = size(endpoints,1);  % edges number
N = 2*(nver+nedge);
M = nver+nedge;
b = zeros(N+M,1);


% -------------------------------------------------------------------------
% Matrices Computation (FEM)
% -------------------------------------------------------------------------

[A, B, C] = RM_fem (xv, yv, vertices, edges, endpoints, ...
                    fdq, mu, lambda, t);


% -------------------------------------------------------------------------
% Matrices Computation (PSRI)
% -------------------------------------------------------------------------

%alpha = max(t^(-2),mu);
ni = 1 / (2*mu + 2*lambda);
alpha = (h^-2)/(5*(1-ni));

if (strcmp(psri,'yes'))
    [A_p, B_p, C_p] = RM_psri (xv, yv, vertices, edges, endpoints, ...
                               fdq, mu, lambda, t, alpha);
    %[A_p, B_p, C_p] = RM_psri_B (xv, yv, vertices, edges, endpoints, ...
    %                             fdq, mu, lambda, t, alpha);
                            
    A = A + A_p;
    B = B + B_p;
    C = C + C_p;
    clear A_p
    clear B_p
    clear C_p
end


% -------------------------------------------------------------------------
% Load Term Computation
% We recall that
% b = [ 0 g]
% where 0 is an N-dimensional array 
% while g is an M-dimensional array 
% -------------------------------------------------------------------------
b(N+1:end) = P2load(xv,yv,vertices,edges,boundary,boundedges,endpoints);


% -------------------------------------------------------------------------
% Border Conditions
% -------------------------------------------------------------------------
internV = setdiff([1:1:nver],boundary); % Indices of internal vertices
internE = setdiff([1:1:nedge],boundedges); % Indices of internal edges
NL = [internV,nver+internE]; % internal dofg (first component)
NL2 = (nver+nedge) + NL;  % internal dofg (second component)
NL3 = (nver+nedge) + NL2;  % internal dofg (third scalar solution)
NL_vect = [NL,NL2]; % internal dofg (first + second component)
NL_load = [NL,NL2,NL3]; % internal dofg (first + second component)

% We now extract the submatrices corrisponding to the internal dofgs
Ah = A(NL_vect,NL_vect);
Bh = B(NL_vect,NL);
BhT = Bh';
Ch = C(NL,NL);
fh = b(NL_load);

% disp(['sizeA = ', num2str(size(A))]);
% disp(['sizeAh = ', num2str(size(Ah))]);
% disp(['sizB = ', num2str(size(B))]);
% disp(['sizeBh = ', num2str(size(Bh))]);
% disp(['sizeb = ', num2str(size(b))]);
% disp(['sizefh = ', num2str(size(fh))]);

% CLear memory from the full matrices that we don't need anymore
clear A
clear B
clear C
clear b


% -------------------------------------------------------------------------
% We now build the final matrix
%         | Ah  |  Bh |
%   Kh =  |-----------|    
%         | BhT |  Ch |
% -------------------------------------------------------------------------
Kh = [Ah Bh; BhT Ch];
uh = zeros(N+M,1);

% disp(['sizeKh = ', num2str(size(Kh))]);
% disp(['sizefh = ', num2str(size(fh))]);

uh(NL_load) = Kh\fh;

% -------------------------------------------------------------------------
% L2 and H1 error computations
% -------------------------------------------------------------------------
[errL2, errH1] = err(uh, fdq, xv, yv, vertices, edges, endpoints, out);

end % end function
