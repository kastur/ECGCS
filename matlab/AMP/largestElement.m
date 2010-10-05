function y = largestElement(x, n)
% LARGESTELEMENT returns nth largest element of x.

% Verify the length of the input
N = length(x);
if(n > N)
    n = N;
    warning('Requested index n is too large, returning smallest element');
elseif(n < 1)
    n = 1;
    warning('Requested index n is too small, returning largest element');

end

% Return the element
t = sort(x, 'descend');
y = t(n);
