function [ z ] = arr_sum( arr )
%Return the sum of all the elements of the array arr

n = length(arr);
z = 0;

for i=1:n
    z = z + arr(i);
end


end

