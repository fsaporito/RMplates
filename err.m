function [errL2, errH1] = err (...
                               uh, ...
                               fdq, ...
                               xv, ...
                               yv, ...
                               vertices, ...
                               edges, ...
                               endpoints, ...
                               out)


% ----------------------------------------------------------------------
% Variables Definition
% ----------------------------------------------------------------------
nver  = length(xv);         % vertices number
nedge = size(endpoints,1);  % edges number
nele = size(vertices,1);    % elements number
N = 2*(nver+nedge);
M = nver+nedge;
theta1_h = uh(1:M);     % theta, first component (size = N/2 = M)
theta2_h = uh(M+1:N); % theta, second component (size = N/2 = M)
w_h = uh(N+1:end);  % w (size = M)
errL2sq = 0;    % L2 error
errH1sq = 0;    % H1 error
uL2 = 0;        % L2 norm exact solution (relative error)
uH1 = 0;        % H1 norm exact solution (relative error)


% ----------------------------------------------------------------------
% Quadrature Formula
% -----------------------------------------------------------------------

% (xhq,yhq) = quadrature node on reference element T_hat
% whq = relative quadrature nodes weight
[xhq, yhq, whq] = quadratura(fdq); % Quadrature Points Computation
Nq = length(xhq); % Number of quadrature points

% Basis functions 
% computed at the quadrature points of the reference element
phihq = zeros(6,Nq);
for i=1:6 % Computation Loop
    for q=1:Nq
        phihq(i,q) = phih2(i,xhq(q),yhq(q));
    end % end for q
end % end for i

% Basis functions gradient 
% computed at the quadrature points of the reference element
gphihqx = zeros(6,Nq);
gphihqy = zeros(6,Nq);
for i=1:6 % Computation Loop
    for q=1:Nq
        [gx,gy] = gradhphih2(i,xhq(q),yhq(q));
        gphihqx(i,q) = gx;
        gphihqy(i,q) = gy;
    end % end for q
end % end for i

% Computation Loop
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
   
    % Solutions in DOFG
    thetaE1 = theta1_h(dofg);   
    thetaE2 = theta2_h(dofg);   
    wE = w_h(dofg);
    
    % TMP Error Variables
    sq  = 0;
    sqH = 0;
    uL = 0;
    uH  = 0;
    
    for q=1:Nq
        
        % Tmp variables for the internal sum
        tmp_t1 = 0; % theta_1 against basis sum
        tmp_t2 = 0; % theta_2 against basis sum
        tmp_w  = 0; % w against basis sum
        tmpH_t1 = [0;0]; % theta_1 against basis gradient sum
        tmpH_t2 = [0;0]; % theta_2 against basis gradient sum
        tmpH_w = [0;0]; % w against basis gradient sum
        
        for i=1:6
            tmp_t1 = tmp_t1 + thetaE1(i)*phihq(i,q);
            tmp_t2 = tmp_t2 + thetaE2(i)*phihq(i,q);
            tmp_w = tmp_w + wE(i)*phihq(i,q);
            tmpH_t1 = tmpH_t1 + thetaE1(i)*(JFIT*[gphihqx(i,q); ...
                                                  gphihqy(i,q)]); 
            tmpH_t2 = tmpH_t2 + thetaE2(i)*(JFIT*[gphihqx(i,q); ...
                                                  gphihqy(i,q)]);
            tmpH_w = tmpH_w + wE(i)*(JFIT*[gphihqx(i,q); ...
                                           gphihqy(i,q)]);
        end % End for i
        
        pos = JF*[xhq(q);yhq(q)]+[x1;y1];
        
        xq = pos(1);
        yq = pos(2);

        % Squared L2 Error
        sq = sq + (theta1(xq,yq) - tmp_t1)^2*whq(q) + ... 
                  (theta2(xq,yq) - tmp_t2)^2*whq(q) + ...
                  (w(xq,yq) - tmp_w)^2*whq(q); 
              
        % Squared H1 Error
        sqH = sqH + norm([theta1_x(xq,yq); 
                          theta1_y(xq,yq)] - tmpH_t1)^2*whq(q) +  ...
                     norm([theta2_x(xq,yq); 
                           theta2_y(xq,yq)] - tmpH_t2)^2*whq(q) +  ...
                     norm([w_x(xq,yq); 
                           w_y(xq,yq)] - tmpH_w)^2*whq(q);
        
        % Exact solution L2 norm (used for relative error)
        uL = uL + (theta1(xq,yq)^2 + ...
                   theta2(xq,yq)^2 + ...
                   w(xq,yq)^2)*whq(q);
               
        % Exact solution H1 norm (used for relative error)      
        uH = uH + (norm([theta1_x(xq,yq); 
                         theta1_y(xq,yq)])^2 + ...
                   norm([theta2_x(xq,yq); 
                         theta2_y(xq,yq)])^2 + ...
                   norm([w_x(xq,yq); 
                         w_y(xq,yq)])^2)*whq(q);       
    end % End for q
    
    % Area Product
    sq = sq*2*area;
    sqH = sqH*2*area;
    uL = uL*2*area;
    uH = uH*2*area;
    
    % Global Error Update
    errL2sq = errL2sq + sq;
    errH1sq = errH1sq + sqH;
    uL2 = uL2 + uL;
    uH1 = uH1 + uH;
    
end % End for iele

% Relative Errors
errL2 = sqrt(errL2sq/uL2);  % Relative L2 Error
errH1 = sqrt(errH1sq/uH1);  % Relative H1 Error

if (strcmp(out,'yes'))
    disp(['      L2 Error (Relative): ' num2str(errL2)]);
    disp(['      H1 Error (Relative): ' num2str(errH1)]);
end % end if

end % end function
