function w_x = w_x(x,y)
%W_X
%    W_X = W_X(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    28-Oct-2017 12:18:56

t2 = y.^2;
t3 = x-1.0;
t4 = t3.^2;
t5 = y-1.0;
t6 = t5.^2;
t7 = x.^2;
t8 = t2-1.0;
t9 = t8.^2;
t10 = t2.^2;
t11 = t4.^2;
t12 = t2.*5.0;
t17 = y.*5.0;
t13 = t12-t17+1.0;
t14 = t7.^2;
t15 = t7.*5.0;
t16 = t15-x.*5.0+1.0;
w_x = t2.*t3.*t4.*t5.*t6.*t7.*y+t2.*t4.*t5.*t6.*t7.*x.*y-t5.*t8.*t9.*t10.*t11.*t13.*t14.*(x.*1.0e1-5.0).*(1.0./1.875e4)-t3.*t4.*t5.*t8.*t9.*t10.*t13.*t14.*t16.*(2.0./9.375e3)-t5.*t7.*t8.*t9.*t10.*t11.*t13.*t16.*x.*(2.0./9.375e3);