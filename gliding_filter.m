function [filtered_A] = gliding_filter(A, varargin);
%[filtered_A] = gliding_filter(A, varargin);
%gliding filter takes 1/3 2/3 1/3 of three adjacent points and puts it in
%the midle.
%Input:
%   A is a matrix and each row will get filtered.
%   varargin indicates how many times the filtering is applied.
%
%Output: 
%   filtered data of matrix A.

if nargin == 1
    namx = 1
else
    nmax = varargin{1};
end

for n = 1:nmax
    d1 = A(:,1);
    d2 = A(:,end);
    A = ([A d2 d2] + [d1 A d2] + [d1 A d2] + [d1 d1 A]);
    A = A(:,2:end-1)/4;
end

filtered_A = A;