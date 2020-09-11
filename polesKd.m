function [poles,residues] = polesKd(funKd,s,r,c,m,n,N,tol)
%POLESKD  Poles of dynamic stiffness matrix.
%   poles = polesKd(funKd,s) returns the vector with the poles of the dynamic
%   stiffness matrix Kd. The matrix-valued function funKd is of dimension s.
%
%   poles = polesKd(funKd,s,r,c) returns the poles in the disk with radius r
%   and center c.
%
%   [poles,residues] = polesKd(...) also returns the vector residuals with the
%   corresponding residues.
%
%   Author:  Roel Van Beeumen
%   Version: 17-09-2015

if nargin < 3, r = 1; end
if nargin < 4, c = 0; end
if nargin < 5, m = 64; end
if nargin < 6, n = 64; end
if nargin < 7, N = 511; end
if nargin < 8, tol = 1e-14; end

% compute poles
z = r*exp(2i*pi*(0:N)/(N+1)) + c;
e = ones(s,1)/sqrt(s);
f = cellfun(@(x) e'*x*e,funKd(z));
[~,~,~,~,~,poles,~] = ratdisk(f,m,n,N,tol);
pol = r*poles + c;

% estimate residues
t = max(tol,1e-7);
residues = 0;
% residues = cellfun(@(x,y) t*((x-y)/2),funKd(pol+t),funKd(pol-t),...
%     'UniformOutput',false);

end % polesKd.m
