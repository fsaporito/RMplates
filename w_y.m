function w_y = w_y(x,y)
%W_Y
%    W_Y = W_Y(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    10-Nov-2017 15:18:29

t2 = x.^2;
t3 = x-1.0;
t4 = t3.^2;
t5 = y.^2;
t6 = y-1.0;
t7 = t6.^2;
t8 = t5.*5.0;
t9 = t8-y.*5.0+1.0;
t10 = t2.*5.0;
t11 = t10-x.*5.0+1.0;
w_y = t2.*t3.*t4.*t6.*t9.*x.*(-1.0./1.75e4)-t3.*t5.*t6.*t7.*t11.*x.*(3.0./1.75e4)-t2.*t3.*t4.*t9.*x.*y.*(1.0./1.75e4)-t3.*t5.*t7.*t11.*x.*y.*(3.0./1.75e4)-t2.*t3.*t4.*t6.*x.*y.*(y.*1.0e1-5.0).*(1.0./1.75e4)+t2.*t3.*t4.*t5.*t6.*t7.*x+t2.*t3.*t4.*t5.*t7.*x.*y;
