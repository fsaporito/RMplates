function z = contr(A,B)
% Return the contraction product between matrices A and B

z = 0;

if (size(A) ~= size(B)) 
  disp('Matrix A:')
  A
  disp('Matrix B:')
  B
  size(A)
  size(B)
  error(message('MATLAB:contr:InputSizeMismatch'));
else
  [N M] = size(A);
  for i=1:N
    for j=1:M
      z = z + A(i,j)*B(i,j);
    end % End For j
  end % End For i

end % End IF

end % End Function