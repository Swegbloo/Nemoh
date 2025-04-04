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
%
%--------------------------------------------------------------------------------------
% This m-file for running NEMOH in MATLAB
% It starts with refining a coarse mesh/aximesh
% It produces Mesh.cal, Mesh.dat, Nemoh.cal
% After meshing done, user suggested to check or adjust the produced
% Nemoh.cal if needed
% After that software runs for NEMOH1 and NEMOH2 (QTF)
%
% User has to specify the output folder in outputdir input below
%--------------------------------------------------------------------------------------
function diff_energy = diffraction_energy(x, w, M)


% x is a vector containing [height, draft, radius]
draft = x(1);
radius_in = x(2);
radius_out = x(3);

% % Define the base path
% basePath = 'C:\Users\swago\OneDrive\Documents\GitHub\Nemoh\matlabRoutines';
% 
% % Add the base path and all its subfolders to the MATLAB path
% addpath(genpath(basePath));


% Dynamically add all subfolders to the MATLAB path
mainFolderPath = 'C:\Users\swago\OneDrive\Documents\GitHub\Nemoh\matlabRoutines';
addpath(genpath(mainFolderPath));

cd('C:\Users\swago\OneDrive\Documents\GitHub\Nemoh\matlabRoutines');

testcase = 4;
ID_PLOT_RESULTS = 1;
ID_QTF = 0; % Flag to
%activate QTF computation (0 or 1)

outputdir = ['.' filesep];
% Update this output files location

% Check that Nemoh is available
assert(FindingNemoh(ID_QTF, true))

switch testcase
    case 1
        %============ MESH WITH Mesh.m==================%
        % RECTANGULAR BOX
        projdir=[outputdir,filesep,'RectangularLiu17'];
        ID_HydrosCal=0;     % A switch,1 computes hydrostatics, inertia, kH.
        % if RAO will be computed, see input in
        % Nemoh.cal, set ID_Hydroscal=1
        
        L = 10; % Length x axis
        W = 20; % Width y axis
        Draft= 5;
        % left face
        X(1,1,:,:)=[-L/2,0,-Draft;
            -L/2,-W/2,-Draft;
            -L/2,-W/2,0;
            -L/2,0,0] ;
        % Bottom face
        X(1,2,:,:)=[-L/2,-W/2,-Draft;
            -L/2,0,-Draft
            L/2,0,-Draft;
            L/2,-W/2,-Draft;];
        % right face
        X(1,3,:,:)=[L/2,0,-Draft;
            L/2,0,0;
            L/2,-W/2,0;
            L/2,-W/2,-Draft] ;
        % front face
        X(1,4,:,:)=[-L/2,-W/2,-Draft;
            L/2,-W/2,-Draft;
            L/2,-W/2,0;
            -L/2,-W/2,0];
        
        ncoarse=length(X(1,:,1,1));   % number of coarse panels
        nBodies=length(X(:,1,1,1));   % 1 body
        tX=0;           % 0 no translation applied to the Mesh
        xyzCoG=[0,0,0]; % position of gravity centre
        nfobj=600;      % target number of panels for Aquaplus mesh
        wmax=3;
        nbfreq=60; % number frequency
        w= linspace(wmax/nbfreq,wmax,nbfreq)'; % min w [rad/s], max w [rad/s], Nw
        wavedir=0;% angle of the incident waves
        depth=600; % water depth (m)
        QTFInput=[ID_QTF,2];%[Switch,Contrib]
        Mesh(nBodies,ncoarse,X,tX,xyzCoG,nfobj,depth,w,wavedir,QTFInput,projdir);
        
    case 2
        %============ MESH WITH Mesh.m==================%
        % RECTANGULAR BOX with lid panels for Ir. freq. removal
        projdir=[outputdir,filesep,'RectangularLiu17_IRR'];
        ID_HydrosCal=1;     % A switch,1 computes hydrostatics, inertia, kH.
        % if RAO will be computed, see input in
        % Nemoh.cal, set ID_Hydroscal=1
        L = 10; % Length x axis
        W = 20; % Width y axis
        Draft= 5;
        % left face
        X(1,1,:,:)=[-L/2,0,-Draft;
            -L/2,-W/2,-Draft;
            -L/2,-W/2,0;
            -L/2,0,0] ;
        % Bottom face
        X(1,2,:,:)=[-L/2,-W/2,-Draft;
            -L/2,0,-Draft;
            L/2,0,-Draft;
            L/2,-W/2,-Draft];
        % right face
        X(1,3,:,:)=[L/2,0,-Draft;
            L/2,0,0;
            L/2,-W/2,0;
            L/2,-W/2,-Draft] ;
        % front face
        X(1,4,:,:)=[-L/2,-W/2,-Draft;
            L/2,-W/2,-Draft;
            L/2,-W/2,0;
            -L/2,-W/2,0];
        % TOP face
        X(1,5,:,:)=[-L/2,-W/2,0;
            L/2,-W/2,0;
            L/2,0,0;
            -L/2,0,0;];
        ncoarse=length(X(1,:,1,1));   % number of coarse panels
        nBodies=length(X(:,1,1,1));   % 1 body
        tX=0;           % 0 no translation applied to the Mesh
        xyzCoG=[0,0,0]; % position of gravity centre
        nfobj=600;      % target number of panels for Aquaplus mesh
        w= linspace(0.05,3,60)'; % min w [rad/s], max w [rad/s], Nw
        wavedir=0;% angle of the incident waves
        depth=600; % water depth (m)
        QTFInput=[ID_QTF,2];%[Switch,Contrib]
        Mesh(nBodies,ncoarse,X,tX,xyzCoG,nfobj,depth,w,wavedir,QTFInput,projdir);
        
    case 3
         %============ MESH WITH axiMesh.m==================%
         %Hemisphere
        projdir=[outputdir,filesep,'hemisphereR10InfD'];
        ID_HydrosCal=0; % A switch,1 computes hydrostatics, inertia, kH.
        % if RAO will be computed, see input in
        % Nemoh.cal, set ID_Hydroscal=1
        
        nang=50;
        npanelt=600;
        
        n1=10; % number of points for the description of the
        n=n1+1;%
        Rayon=10;
        for i=1:n+1
            r(i)=Rayon*cos((n1-i+1)*pi/2/n1-pi/2);
            z(i)=Rayon*sin((n1-i+1)*pi/2/n1-pi/2);
        end
        
        %------------Define calculation options-------------
        g=9.811;
        w= linspace(0.05,3,60)'; % min w [rad/s], max w [rad/s], Nw
        wavedir=0;               % angle of the incident waves
        depth=600;               % water depth (m)
        xyzCoG=[0,0,0];          % position of gravity centre
        QTFInput=[ID_QTF,2];     %[Switch,Contrib]
        
        axiMesh(r,z,n,nang,npanelt,xyzCoG,depth,w,wavedir,QTFInput,projdir);% Call the function axiMesh.m

    case 4
         %============ MESH WITH axiMesh.m==================%
         %Cylinder
        projdir=[outputdir,filesep,'cylinder'];
        ID_HydrosCal=1; % A switch,1 computes hydrostatics, inertia, kH.
        % if RAO will be computed, see input in
        % Nemoh.cal, set ID_Hydroscal=1
        
        QTFInput=[ID_QTF,2]; 

        % Verify if Nemoh.cal exists
    if ~isfile([projdir, filesep, 'Nemoh.cal'])
        error('Nemoh.cal was not generated in the directory: %s', projdir);
    end
end

%% meshing
generate_cylinder_mesh(draft, radius_in, radius_out, M); % mesh making and directly edit geometry file
    
    
disp('Meshing done!')
disp('Check and adjust the produced Nemoh.cal in the project directory!')
disp('Press enter for the NEMOH computation')
%pause;
%-------Launch Calculation------------
%[Idw,w,A,B,Fe]=Nemoh(projdir,ID_HydrosCal,ID_QTF); % Call the function Nemoh.m


%% solving for diffraction energy
[nw] = Nemoh(projdir,ID_HydrosCal,ID_QTF,w); % Call the function Nemoh.m
kochmag = zeros(1,nw);
kochang = zeros(1,nw);
diff_energy = zeros(nw,1);
%% --- Computes QTFs --------------------
if ID_QTF==1, NemohQTF(projdir);end % Call the function NemohQTF.m

%% ---------------------------------------
% if ID_PLOT_RESULTS==1 %Plot NEMOH1
%     if Idw==1
%         xlab='\omega [rad/s]';
%     elseif Idw==2
%         xlab='f [Hz]';
%     elseif Idw==3
%         xlab='T [s]';
%     end
% 
%     figure
%     subplot(3,2,1)
%     a(1,:)=A(1,1,:);
%     b(1,:)=B(1,1,:);
%     plot(w,a,'b-+',w,b,'r-+')
%     grid ON
%     ylabel('Coeffs' )
%     legend('Added mass (A)','Damping (B)')
%     title('Surge')
%     plot_properties;
% 
%     subplot(3,2,2)
%     plot(w,abs(Fe(:,1)),'g-+')
%     grid ON
%     ylabel('F (N)' )
%     title('Surge')
%     plot_properties;
% 
%     subplot(3,2,3)
%     a(1,:)=A(3,3,:);
%     b(1,:)=B(3,3,:);
%     plot(w,a,'b-+',w,b,'r-+')
%     grid ON
%     ylabel('Coeffs' )
%     legend('Added mass (A)','Damping (B)')
%     title('Heave')
%     plot_properties;
% 
%     subplot(3,2,4)
%     plot(w,abs(Fe(:,3)),'g-+')
%     grid ON
%     ylabel('F(N)' )
%     title('Heave')
%     plot_properties;
% 
%     subplot(3,2,5)
%     a(1,:)=A(5,5,:);
%     b(1,:)=B(5,5,:);
%     plot(w,a,'b-+',w,b,'r-+')
%     xlabel(xlab)
%     ylabel('Coeffs' )
%     grid ON
%     legend('Added mass (A)','Damping (B)')
%     title('Pitch')
%     plot_properties;
% 
%     subplot(3,2,6)
%     plot(w,abs(Fe(:,5)),'g-+')
%     xlabel(xlab)
%     grid ON
%     ylabel('F(N)' )
%     title('Pitch')
%     plot_properties;
% 
%     if ID_QTF==1
%         qtftype='M';
%         Idwvar=0;%0=w rad/s,1=f 1/s,2=T s; %freq type
%         SwitchBiDir=0;
%         betaID=[0 0]; %[beta1 beta2]
%         shiftdw1=-1;
%         ShowLegend=1;
%         Idaxisij=0;
%         NWdatNEM=length(w);%specify number of freq as in the QTF output
%         DOFdatNem=[1 1 1 1 1 1];% the QTF data availability in each mode, 1 means data exist
%         % if for index 1 to 6, all the values are 1, it
%         % means data available for all DOF
%         DOFsurge=1;DOFheave=3;DOFpitch=5;
%         NbetaData=1;
% 
%         Id_fixed_CB_Axis=0; %switch to fixed colorbar limit, 
%                             % if 1 the minmaxQ values have to be specified
%         minmaxQR_surge=[-20 20]; minmaxQI_surge=[-20 20];minmaxQMod_surge=[0 20];
%         minmaxQR_heave=[-10 10]; minmaxQI_heave=[-10 10];minmaxQMod_heave=[0 10];
%         minmaxQR_pitch=[-50 50];minmaxQI_pitch=[-50 50];minmaxQMod_pitch=[0 50];
% 
%         plot_QTF_NEMOH;
%     end
% end
%fid = fopen([projdir, filesep, 'results', filesep, 'Kochin.00001.dat'], 'r');

%% extract energies

% Define the range of numbers for the files
startNum = 1; % Start from 1 (for Kochin.00001.dat)
endNum = nw;   % 'n' is the total number of files (<100)

% Loop through each file
for i = startNum:endNum
    % Construct the file name dynamically
    fileName = sprintf('Kochin.%05d.dat', i); % '%05d' ensures zero padding to 5 digits
    filePath = [projdir, filesep, 'results', filesep, fileName]; % Full path

    % Open the file
    fid = fopen(filePath, 'r');
    
    % % Check if the file was opened successfully
    % if fid == -1
    %     fprintf('Could not open %s. Skipping...\n', fileName);
    %     continue; % Skip to the next file
    % end

    % Process the file (example: read the first line)
    ligne = fgetl(fid);
    %disp(['First line of ', fileName, ': ', ligne]);

    % Close the file
    fclose(fid);
    % Open file and read line if from file
    % fid = fopen('yourfile.dat', 'r'); % Uncomment if reading from file
    % ligne = fgetl(fid);               % Read line from file
    
    % Split the line into numbers using sscanf
    numbers = sscanf(ligne, '%f'); % Parse the line into an array of floats
    
    % Extract the second number
    kochmag(1,i) = numbers(2);

    % Extract the third number
    kochang(1,i) = numbers(3);
    
    % Display the result
    %disp(['The third number is: ', num2str(thirdNumber)]);
    
    % fclose(fid); % Close the file when done
    
    % [minE,k_0] = praxis(0.01,2,3,0,1,@);
    diff_energy(i,1) = -cos(kochang(1,i))*2*kochmag(1,i);
end


