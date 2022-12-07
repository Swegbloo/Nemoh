%--------------------------------------------------------------------------------------
%
%    Copyright (C) 2022 - LHEEA Lab., Ecole Centrale de Nantes, UMR CNRS 6598
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   Contributors list:
%   - R. Kurnia
%--------------------------------------------------------------------------------------
%
% --> function []=NemohQTF(bindir,projdir)
%
% Purpose: NEMOH Matlab wrapper for computing QTFs
%
% Inputs :
% - bindir : Binary directory
% - projdir: Project directory
%
%
function NemohQTF(bindir,projdir)

system(['mkdir ',projdir,filesep,'QTFPreprocOut']);
system(['mkdir ',projdir,filesep,'results',filesep,'QTF']);

% Calcul des coefficients hydrodynamiques
l = isunix;
if l == 1
    fprintf('\n------ Starting QTF preproc ----------- \n');
    system([bindir,filesep,'QTFpreProc ',projdir]);
    fprintf('------ Computing QTFs ------------- \n');
    system([bindir,filesep,'QTFsolver ',projdir]);
    fprintf('------ Postprocessing results --- \n');
    system([bindir,filesep,'QTFpostProc ',projdir]);
else
    fprintf('\n------ Starting QTF preproc ----------- \n');
    system([bindir,filesep,'QTFpreProc.exe ',projdir]);
    fprintf('------ Computing QTFs ------------- \n');
    system([bindir,filesep,'QTFsolver.exe ',projdir]);
    fprintf('------ Postprocessing results --- \n');
    system([bindir,filesep,'QTFpostProc.exe ',projdir]);
end

end
