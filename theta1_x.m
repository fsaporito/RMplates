function theta1x = theta1_x(x,y)
%THETA1_X
%    THETA1X = THETA1_X(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 5.8.
%    10-Nov-2017 15:18:28

t2 = y.^2;
t3 = x-1.0;
t4 = y-1.0;
t5 = t4.^2;
t6 = t3.^2;
t7 = x.^2;
t8 = x.*2.0;
t9 = t8-1.0;
theta1x = t2.*t4.*t5.*t6.*t7.*y.*2.0+t2.*t4.*t5.*t6.*t9.*x.*y.*2.0+t2.*t4.*t5.*t7.*t9.*y.*(t8-2.0);
