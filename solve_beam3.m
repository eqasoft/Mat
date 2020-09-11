% /*****************************************************************************
%   Copyright (c) 2020, EQASoft IN.
%   All rights reserved.
% 
%   Redistribution and use in source and binary forms, with or without
%   modification, are permitted provided that the following conditions are met:
% 
%     * Redistributions of source code must retain the above copyright notice,
%       this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of EQASoft IN. nor the names of its contributors
%       may be used to endorse or promote products derived from this software
%       without specific prior written permission.
% 
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
%   THE POSSIBILITY OF SUCH DAMAGE.
% 
% ******************************************************************************
% * Contents:  This file contains the 'main' function. Program execution begins
%              and ends there.
% 
% * Developed by :
%             Avisek Mukherjee (amukherjee@eqasoft.com)
% 
% * Generated July 2015
% *****************************************************************************/


tic;
% solve_beam
clear all;
%close all;
clc;
%% Structural Properties
nElem = 20;
Lx = 10;
Material{1}.EI = ones(1,nElem)*15.625E7;             % N-m^2 
Material{1}.EA = ones(1,nElem)*7.5E9;                % N 
Material{1}.nu = 0.3;
Material{1}.m = ones(1,nElem)*625;                   % kg/m
Material{1}.xi = 0.0;                                % 0% of material damping

u1 = 0;
u2 = 0;


%% 1D KL decomposition
klorder =10;
Lc = 1;
                                              % Stochastic domain length
%E_mean = 30*10^9;                                   % Mean value of Young's modulus (stochastic process E)
%E_std = 0.3*E_mean;                                 % Standard deviation of Young's modulus


% 1D Beam Element
Ex = linspace(0,Lx,nElem+1)';

% mesh
Nodes = [(1:(nElem+1))' Ex zeros(nElem+1,1) zeros(nElem+1,1)];


cp = size(Nodes,1);
C = zeros(cp);
for j=1:cp
   for k=1:cp
      C(j,k) = exp(-sqrt(sum(((Nodes(j,2:4)-Nodes(k,2:4)).^2)))/Lc);  
   end
end 

[phi,lambda] = eig(C);
lambda = diag(lambda);
[lambda,ilambda] = sort(lambda,1,'descend');
lambda_norm = diag(lambda)/cp;
phi = phi(:,ilambda);

lambda = lambda(1:klorder,1);
phi = phi(:,1:klorder);





%% Pick KL Mode Values

dx = Lx/nElem;
%sampleP = zeros(nElem,klorder);
for iKL = 1:klorder
    sampleP = interp1q(Ex,phi(:,iKL),Ex(1:(end-1))+dx/2);
    
    for iElem = 1:nElem
        Material{iElem+1}.EI = sampleP*15.625E7*0.3;                % N-m^2 
        Material{iElem+1}.EA = sampleP*7.5E9*0.3;                % N 
        Material{iElem+1}.nu = 0.3;
        Material{iElem+1}.m  = sampleP*625*0.3;                  % kg/m
        Material{iElem+1}.xi = 0.0;                % 5% of material damping
    end
    
end




%%
% Standard normal random variable
rng('default') % For reproducibility
%Kd = klCompose(f,nElem,D,Lx,klorder);
Nsample = 1000000; 
KL_xi_sample = normrnd(0,1,[klorder,Nsample]);
%%

numIterations = length(100);
ppm = ParforProgressbar(numIterations, 'parpool', {'local', 16});
%% Function Handler for spectral element
parfor iSample = 2201:2300
    KL_xi = KL_xi_sample(:,iSample);
    funKd = @klCompose;
    fn   = 50;        % No. of samples in frequency domain
    fmin = 0.001;   % Minimum frequency (Hz)
    fmax = 400;    % Maximum frequency (Hz)

    % interpolation points on circle
    radius = (fmax - fmin)/2;
    center = fmin + radius;
    z = exp(2i*pi*(0:fn-1)/(fn));
    f = radius*z + center;

    % samples of Kd
    Kd = funKd(f,nElem,Material,Lx,klorder,KL_xi,lambda);
    %Kd = funKd(f,nElem,Material,Lx,KL_xi,lambda);
    s = size(Kd{1},1);


    % pole handling
    e = ones(s,1)/sqrt(s);
    ff = cellfun(@(x) e'*x*e,funKd(f,nElem,Material,Lx,klorder,KL_xi,lambda));
    [~,a,b,~,~,~,~] = ratdisk(ff,floor((fn-1)/2),floor((fn-1)/2),(fn-1),1e-14);
    poles = roots(b(end:-1:1));

    % remove poles
    P = Kd;
    for i = 1:length(P)
        for j = 1:length(poles)
            P{i} = P{i}*(z(i) - poles(j));
        end
        P{i} = P{i};
    end

    % scaling
    maxnrmP = max(cellfun(@(x) norm(x),P));
    for i = 1:length(P)
        P{i} = P{i}/maxnrmP;
    end


     % solve nlep
    [C0,C1] = linlagr(P,z,z);   % linearization
    [V,D] = eig(C0,C1);         % solve glep

    lam = radius*diag(D) + center;
    err = zeros(size(lam));

    for i = 1:length(lam)
        u = V(1:s,i);
        u = u/norm(u);
        Kdlam = funKd(lam(i),nElem,Material,Lx,klorder,KL_xi,lambda);
        err(i) = norm(Kdlam*u);
    end


    % output
    tol = 1e-6;
    sol{iSample} = lam(err<tol);


    % save solve_frame_linearization;
    % plot
%     figure
%     hold off; plot(radius*exp(2i*pi*(0:fn-1)/(fn)) + center,'-k');
%     hold on;
%     plot(lam,'*b');
%     plot(lam(err<tol),'*b');
%     plot(radius*poles + center+1i*1E-10,'or');
%     xlim([fmin fmax])
%     ylim([-fmax/2 fmax/2])

    %sol =   sol(unique(round(sort(real(sol)),4))>=0);
    feq{iSample} = unique(round(sort(real(sol{iSample})),4));
    feq{iSample} = feq{iSample}(feq{iSample}>=0);
    
    
 ppm.increment();
end
delete(ppm)

save result_beam_KL_2201to2300.mat;

for index = 1:size(feq,1)
    size_feq(index) = size(feq{index},1);
end
min(size_feq)

toc;
