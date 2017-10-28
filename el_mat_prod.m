function z = el_mat_prod(A,B)
% Return the element wise product between matrices A and B

z = 0;

if (size(A) ~= size(B)) 
  disp('Matrix A:')
  A
  disp('Matrix B:')
  B
  size(A)
  size(B)
  error(message('MATLAB:mat_el_prod:InputSizeMismatch'));
else
  [N M] = size(A);
  z = zeros(N,M);
  for i=1:N
    for j=1:M
      z(i,j) = A(i,j)*B(i,j);
    end % End For j
  end % End For i

end % End IF

end % End Function