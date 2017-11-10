function [] = pre_comp_sym (m_u,l_u,t) 

syms x y real;

% Lambda
lambda = l_u;

% Mu
mu = m_u;

% Ni
ni = 1 / (2*mu + 2*lambda);

% E
E = mu/(lambda + mu);

% First Component
%theta1 = cos(x)*(x^2 + y^2)
% theta1 = x^4;
theta1 = y^(3) * ( y - 1 )^(3) * x^(2) * ( x - 1 )^(2) * ( 2*x - 1 );

% Second Component
%theta2 = sin(x)*(x^2 + y^2)
% theta2 = y^4;
theta2 = x^(3) * ( x - 1 )^(3) * y^(2) * ( y - 1 )^(2) * ( 2*y - 1 );


% First Derivatives, First Component
theta1x = diff(theta1,x);
theta1y = diff(theta1,y);

% Second Derivatives, First Component
theta1xx = diff(theta1x,x);
theta1xy = diff(theta1x,y);
theta1yx = diff(theta1y,x);
theta1yy = diff(theta1y,y);

% First Derivatives, Second Component
theta2x = diff(theta2,x);
theta2y = diff(theta2,y);

% Second Derivatives, Second Component
theta2xx = diff(theta2x,x);
theta2xy = diff(theta2x,y);
theta2yx = diff(theta2y,x);
theta2yy = diff(theta2y,y);

% Partial Derivative of w in x
f1 = (1/24)*(lambda*theta1xx + ...
             lambda*theta2xy + ...
             mu*theta1xx + ...
             mu*theta2xy + ...
             mu*theta1yy);
    
%w_x = theta1 - (t^2)/(24*mu) * f1;
        
% Partial Derivative of w in  y
f2 = (1/24)*(mu*theta2xx + ...
             mu*theta1yx + ...
             lambda*theta1xy + ...
             lambda*theta2yy + ...
             mu*theta1xy);
    
%w_y = theta2 - (t^2)/mu) * f2;
        
 % w Computation as Potential 
 % w = potential([w_x w_y], [x y]);
 w = (1/3)*x^(3) * (x-1)^(3) * y^(3) * (y-1)^(3) ...
     - (2/5)*(t^(2)/(1-ni)) * ...
     (y^(3) * (y-1)^(3) * x * (x-1) * (5*x^(2) - 5*x + 1) + ...
      x^(3) * (x-1)^(3) * y * (y-1) * (5*y^(2) - 5*y + 1) );
  w_x = diff(w,x);
  w_y = diff(w,y);
    
% Load
% g = -div (t^(-2)*mu*(nabla w - theta))
g1 = -diff( f1, x) + ...
    -diff( f2, y);
g2 = -t^(-2)*mu*(diff( w_x - theta1, x) + ...
    		 diff( w_y - theta2, y));
g3 = (E / ( 12 * ( 1 - ni^(2) ) ) ) * ...
     ( 12*y * ( y - 1 ) * ( 5*x^(2) - 5*x + 1 ) * ...
		 ( 2*y^(2) * ( y - 1 )^(2) + ...
		   x * ( x - 1 ) * ( 5*y^(2) - 5*y + 1 )) + ...
       12*x * (x-1) * ( 5*y^(2) - 5*y + 1 ) * ...
         ( 2*x^(2) * ( x - 1 )^(2) + ... 
           y * ( y - 1 ) * ( 5*x^(2) - 5*x + 1 )) ...
      );


% Writing functions on files
% 'Vars' [x y] forces the input in both  x and y
%matlabFunction(lambda,'File','lambda');
%matlabFunction(mu,'File','mu');
matlabFunction(theta1,'File','theta1','Vars',[x y]);
matlabFunction(theta1x,'File','theta1_x','Vars',[x y]);
matlabFunction(theta1y,'File','theta1_y','Vars',[x y]);
matlabFunction(theta2,'File','theta2','Vars',[x y]);
matlabFunction(theta2x,'File','theta2_x','Vars',[x y]);
matlabFunction(theta2y,'File','theta2_y','Vars',[x y]);
matlabFunction(w,'File','w','Vars',[x y]);
matlabFunction(w_x,'File','w_x','Vars',[x y]);
matlabFunction(w_y,'File','w_y','Vars',[x y]);
matlabFunction(g3,'File','load_plate','Vars',[x y]);
matlabFunction(g1,'File','load_plate1','Vars',[x y]);
matlabFunction(g2,'File','load_plate2','Vars',[x y]);
matlabFunction(g3,'File','load_plate3','Vars',[x y]);


%reset(symengine)