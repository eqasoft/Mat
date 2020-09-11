function [C0,C1] = linlagr(A,x,w)
%LINLAGR  Linearization of Lagrange matrix polynomial.
%   [C0,C1] = LINLAGR(A,x,w) returns the companion pencil (C0,C1) obtained from
%   the linearization of the interpolating (barycentric) Lagrange polynomial. A
%   is a cell array containing the interpolation values in the interpolation
%   points x. The corresponding barycentric weights are given by w.
%
%   References:
%   R. Van Beeumen, W. Michiels, K. Meerbergen, Linearization of Lagrange and
%      Hermite interpolating matrix polynomials. 2013
%
%   Author:  Roel Van Beeumen
%   Version: 12-05-2013

n = length(x) - 1;         % polynomial degree n
s = size(A{1},1);          % polynomial dimension s
t = w(1:end-1)./w(2:end);  % ratios of barycentric weights

%% matrices C0 & C1
C0 = zeros(n*s);
C1 = zeros(n*s);

% first row
blk1 = 1:s;
blkn = (n-1)*s+1:n*s;
for i = 1:n-1
    blki = (i-1)*s+1:i*s;
    C0(blk1,blki) = x(i+1)*A{i};
    C1(blk1,blki) = A{i};
end
C0(blk1,blkn) = x(n+1)*A{n} + x(n)/t(n)*A{n+1};
C1(blk1,blkn) = A{n} + 1/t(n)*A{n+1};

% diagonal and sub-diagonal
for i = 1:n-1
    blki = (i-1)*s+1:i*s;
    % diagonal
    C0(blki+s,blki+s) = -x(i+2)*t(i)*eye(s);
    C1(blki+s,blki+s) = -t(i)*eye(s);
    % sub-diagonal
    C0(blki+s,blki) = x(i)*eye(s);
    C1(blki+s,blki) = eye(s);
end
