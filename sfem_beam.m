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
% `
% ******************************************************************************
% * Contents:  This file contains the 'main' function. Program execution begins
%              and ends there.
% 
% * Developed by :
%             Avisek Mukherjee (amukherjee@eqasoft.com)
% 
% * Generated July 2015
% ************************************************************************/

% Computation of frequency response function of a cantilever beam modeled
% with two spectral beam elements of 3m and 7m long respectively.

% Avisek Mukherjee
% July 2015

%-------------------------------------------------------------------------
function Kd = sfem_beam(f,EI,EA,m,L)

omega=2*pi*f;

% Material properties
% D.EI = 15.625E7;                % N-m^2 
% D.EA = 7.5E9;                % N 
% D.m = 625;                  % kg/m

D.nu = 0.2;
D.xi = 0.0;                % Material damping

D.EI  = EI;
D.EA = EA;                % N 
D.m = m;                  % kg/m

% Element coefficient matrices
A1 = zeros(6);
A1(2,2) = -1; 
A1(5,5) = -1;

A2 = zeros(6);
A2(2,5) = 1; 
A2(5,2) = 1; 

A3 = zeros(6);
A3(3,2) = -1;    
A3(2,3) = -1;  
A3(5,6) = 1;  
A3(6,5) = 1;   

A4 = zeros(6);    
A4(2,6) = -1; 
A4(6,2) = -1; 
A4(3,5) = 1; 
A4(5,3) = 1;    

A5 = zeros(6);    
A5(6,6) = -1; 
A5(3,3) = -1;  

A6 = zeros(6);
A6(3,6) = -1;
A6(6,3) = -1;

A8 = zeros(6);
A8(1,1) = 1; 
A8(4,4) = 1; 

A9 = zeros(6);
A9(1,4) = -1; 
A9(4,1) = -1;



% Element-1 
%L=4; 
for ifreq = 1:1
    
    % Along 'Z' axis : Flexural
    kz = ((omega(ifreq).^2 * D.m )/(D.EI))^(1/4);

    cz = cos(kz*L); 
    sz = sin(kz*L) ;
    Cz = cosh(kz * L);
    Sz = sinh(kz*L);
              
    f1 = (Cz * sz + Sz * cz);   
    f2 = (Sz + sz);    
    f3 = (Sz * sz);             
    f4 = (Cz - cz);    
    f5 = (Cz * sz - Sz * cz);   
    f6 = (Sz - sz);    
    f7 = (Cz * cz - 1);
    
    % Along 'X' axis : Axial
    kx = ((omega(ifreq).^2 * D.m )/(D.EA))^(1/2);
    cx = cos(kx*L); 
    sx = sin(kx*L);

    f8 = cx / sx;               
    f9 = 1/ sx; 

    Kd = (D.EI/f7)*((kz^3 * (A1*f1 + A2*f2)) + (kz^2 * (A3*f3 + A4*f4))...
        +(kz * (A5*f5 + A6*f6))) + (D.EA*kx*(A8*f8 + A9*f9));
    Kd = Kd/1e9;
end

