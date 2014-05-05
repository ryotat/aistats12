function x = randstgs(l, u, x)
% RANDTG Random numbers from standard truncated Gaussian density
%   p(x) = const * exp(-x^2/2) *I(l < x < u) using slice sampling
%
% Usage
%   x = randtg(l, u, x0) returns an array of random numbers chosen 
%     from the truncated Gaussian density. The size of x is the maximum 
%     common size of the parameters.
%
% Input
%   l, u         Truncation, l < x < u

% Copyright 2009 Mikkel N. Schmidt, ms@it.dk, www.mikkelschmidt.dk

sz = accumarray([1:ndims(l),1:ndims(u),1:ndims(x)]',[size(l),size(u),size(x)]',[],@max)';

z = bsxfun(@plus, -.5*x.^2, log(rand(sz)));
s = sqrt(-2*z);
ll = bsxfun(@max, -s, l);
uu = bsxfun(@min, s, u);
x = rand(sz).*(uu-ll)+ll;