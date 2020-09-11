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
% * Generated September, 2019
% *****************************************************************************/


function P = removepoles(Kd,f,poles,residues)



% remove poles
P = Kd;
for i = 1:length(P)
    for j = 1:length(poles)
        P{i} = P{i}*(f(i) - poles(j));
    end
    P{i} = P{i};
end

% scaling
maxnrmP = max(cellfun(@(x) norm(x),P));
for i = 1:length(P)
    P{i} = P{i}/maxnrmP;
end

end % removepoles.m
