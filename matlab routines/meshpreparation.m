clc
clear all
close all

testcase=12; % 1: Rectangular, 2: Cylinder, 3: half submerge-sphere, 4: HELOFOW monopile
            % 5: Hemisphere KIm Yue QTF 6: Cylinder KIm Yue QTF D/R=4 7:Cylinder KIm Yue QTF D/R=1
            % 8: Softwind nang50 npt 800 9:  Softwindnang 50 npt 1200
switch testcase 
%     case 1
%         %============ MESH WITH Mesh.m==================%
%         % RECTANGULAR BOX
%         clc
%         clear all
%         close all
%         L = 10; % Length x axis
%         W = 20; % Width y axis
%         Draft= 5;
%         % left face
%         X(1,1,:,:)=[-L/2,0,-Draft;
%             -L/2,-W/2,-Draft;
%             -L/2,-W/2,0;
%             -L/2,0,0] ;
%         % Bottom face
%         X(1,2,:,:)=[-L/2,-W/2,-Draft;
%             -L/2,0,-Draft
%             L/2,0,-Draft;
%             L/2,-W/2,-Draft;];
%         % right face
%         X(1,3,:,:)=[L/2,0,-Draft;
%             L/2,0,0;
%             L/2,-W/2,0;
%             L/2,-W/2,-Draft] ;
%         % front face
%         X(1,4,:,:)=[-L/2,-W/2,-Draft;
%             L/2,-W/2,-Draft;
%             L/2,-W/2,0;
%             -L/2,-W/2,0];
%         
%         n=4;         % number of panels
%         nBodies=1;   % 1 body
%         tX=0;        % no translation applied to the Mesh
%         CG=[0,0,-2.5]; % position of gravity centre
%         nfobj=350;   % target number of panels for Aquaplus mesh
%         [Mass,Inertia,KH,XB,YB,ZB]=Mesh(nBodies,n,X,tX,CG,nfobj);
%         M(:,:)=Inertia(1,:,:); % Mass Matric of the first body
%         KHyd(:,:)=KH(1,:,:);   % Hydrostatic stiffness of the first body
%         
%         save('Mesh_outputs','KHyd','M')
%         
%          %------------Define calculation options-------------
%         nbfreq=10; % number of calculations = (number of BVP per frequency*number of frequencies)
%         w= linspace(2*pi/20,2*pi/3,nbfreq)'; % Periods from 3s to 20s for waves for instance
%         dir=0;% angle of the incident waves
%         depth=60; % water depth (m)
%         
%         %============ MESH WITH AXIMESH.M=================%
%         %% VERTICAL CYLINDER
%     case 2
%         clc
%         clear all
%         close all
%         
%         n=3; % 3 points are required for describing the shape
%         Radius =5; % Radius of the cylinder
%         Draft=-10; % Height of the submerged part
%         r=[Radius Radius 0]; % r is the first coordinates of the
%         z=[0 Draft Draft];
%         [Mass,Inertia,KHyd,XB,YB,ZB]=axiMesh(r,z,n);% Call the function axiMesh.m
%         M=Inertia; % Mass Matrix
%         
%         save('Mesh_outputs','KHyd','M')
%         %% Half-submerged Spherel
%         
%          %------------Define calculation options-------------
%         nbfreq=10; % number of calculations = (number of BVP per frequency*number of frequencies)
%         w= linspace(2*pi/20,2*pi/3,nbfreq)'; % Periods from 3s to 20s for waves for instance
%         dir=0;% angle of the incident waves
%         depth=60; % water depth (m)
%     case 3
%         %  Don't Forget z(i)> z(i+1) !!
%         clc
%         clear all
%         close all
%         
%         n1=10; % number of points for the description of the
%         n=n1+1; %
%         Rayon=1;
%         for i=1:n+1
%             r(i)=Rayon*cos((n1-i+1)*pi/2/n-pi/2);
%             z(i)=Rayon*sin((n1-i+1)*pi/2/n-pi/2);
%         end
%         [Mass,M,KHyd,XB,YB,ZB]=axiMesh(r,z,n); % M is the Mass Matrix, and KHyd the Hydrostatic stiffness
%         
%         save('Mesh_outputs','KHyd','M') % For reusing the 2 matrices in the RAO calculation
%         
%         %------------Define calculation options-------------
%         nbfreq=10; % number of calculations = (number of BVP per frequency*number of frequencies)
%         w= linspace(2*pi/20,2*pi/3,nbfreq)'; % Periods from 3s to 20s for waves for instance
%         dir=0;% angle of the incident waves
%         depth=60; % water depth (m)
% 
%         %%
%     case 4
%             %% HelloFOW Spar 
%             %% specify n angular 50, npanel 500, z gravity -1.8*Frscaled
%         clc
%         clear all
%         close all
%         
%         n=5; % 3 points are required for describing the shape
%         FrScaled=40;
%         R1 =0.14.*FrScaled;
%         R2 =0.225.*FrScaled; % Radius of the cylinder
%         Draft =-2.25.*FrScaled; % Height of the submerged part
%         r=[R1 R1 R2 R2 0]; % r is the first coordinates of the
%         z=[0 -0.137.*FrScaled -0.337.*FrScaled Draft Draft];
%         [Mass,Inertia,KHyd,XB,YB,ZB]=axiMesh(r,z,n);% Call the function axiMesh.m
%         M=Inertia; % Mass Matrix
%         
%         save('Mesh_outputs','KHyd','M')
%         
%         nbfreq=10; % number of calculations = (number of BVP per frequency*number of frequencies)
%         w= linspace(2*pi/20,2*pi/3,nbfreq)'; % Periods from 3s to 20s for waves for instance
%         dir=0;% angle of the incident waves
%         depth=5*FrScaled; % water depth (m)
%     case 5 %Hemisphere KIM&Yue QTF testcase
%          
%          %  Don't Forget z(i)> z(i+1) !!
%         clc
%         clear all
%         close all
%         
%         n1=21; % number of points for the description of the
%         n=n1+1; %
%         Rayon=10; %radius
%         %zz=linspace(0,Rayon,n+1);
%         for i=1:n+1
%             r(i)=Rayon*cos((n1-i+1)*pi/2/n1-pi/2);
%             z(i)=Rayon*sin((n1-i+1)*pi/2/n1-pi/2);
%         end
%        
%         [Mass,M,KHyd,XB,YB,ZB]=axiMesh(r,z,n); % M is the Mass Matrix, and KHyd the Hydrostatic stiffness
%        % axis equal
%         save('Mesh_outputs','KHyd','M') % For reusing the 2 matrices in the RAO calculation
%         
%         %------------Define calculation options-------------
%         nbfreq=40; % number of calculations = (number of BVP per frequency*number of frequencies)
%         w= linspace(2/nbfreq,2,nbfreq)'; % Periods from 3s to 20s for waves for instance
%         dir=0;% angle of the incident waves
%         depth=30; % water depth (m)
%   case 6 %Cylinder KIM&Yue QTF test
%         clc
%         clear all
%         close all
%         
%         n=3; % 3 points are required for describing the shape
%         Radius =10; % Radius of the cylinder
%         Draft=-10; % Height of the submerged part
%         r=[Radius Radius 0]; % r is the first coordinates of the
%         z=[0 Draft Draft];
%         [Mass,Inertia,KHyd,XB,YB,ZB]=axiMesh(r,z,n);% Call the function axiMesh.m
%         M=Inertia; % Mass Matrix
%         
%         save('Mesh_outputs','KHyd','M')
%         %% Half-submerged Spherel
%         
%          %------------Define calculation options-------------
%         nbfreq=40; % number of calculations = (number of BVP per frequency*number of frequencies)
%         w= linspace(2/nbfreq,2,nbfreq)'; % Periods from 3s to 20s for waves for instance
%         dir=0;% angle of the incident waves
%         depth=40; % water depth (m)
%         
%     case 7 %Cylinder KIM&Yue QTF test D/R=1
%         clc
%         clear all
%         close all
%         
%         n=3; % 3 points are required for describing the shape
%         Radius =10; % Radius of the cylinder
%         Draft=-10; % Height of the submerged part
%         r=[Radius Radius 0]; % r is the first coordinates of the
%         z=[0 Draft Draft];
%         [Mass,Inertia,KHyd,XB,YB,ZB]=axiMesh(r,z,n);% Call the function axiMesh.m
%         M=Inertia; % Mass Matrix
%         
%         save('Mesh_outputs','KHyd','M')
%         %% Half-submerged Spherel
%         
%          %------------Define calculation options-------------
%         nbfreq=40; % number of calculations = (number of BVP per frequency*number of frequencies)
%         w= linspace(2/nbfreq,2,nbfreq)'; % Periods from 3s to 20s for waves for instance
%         dir=0;% angle of the incident waves
%         depth=10; % water depth (m)     
        
    
        case 8
            %% SOFTWIND PLATFORM
            %% specify n angular 50, npanel 500, z gravity -1.789*Frscaled
        clc
        clear all
        close all
        
        dirname='SOFTWINDNang50Npt800';
        FrScaled=40;
        R1 =0.14.*FrScaled;
        R2 =0.225.*FrScaled; % Radius of the cylinder
        R121=0.1683*FrScaled;
        R122=0.1966*FrScaled;
        
        Draft =-2.285.*FrScaled; % Height of the submerged part
        
        r=[R1 R1 R121 R122 R2 R2 3*R2/4 2*R2/4 1*R2/4 0]; % r is the first coordinates of the
        z=[0 -0.135.*FrScaled -0.2016*FrScaled -0.2682*FrScaled -0.335.*FrScaled Draft Draft Draft Draft Draft];
        n=length(r);
        zCoG=-1.789*FrScaled;
        nang=50;
        npanelt=800;
        nbfreq=60; % number of calculations = (number of BVP per frequency*number of frequencies)
        w= linspace(0.05,3,nbfreq)'; % Periods from 3s to 20s for waves for instance
        dir=0;% angle of the incident waves
        depth=5*FrScaled; % water depth (m)
        QTFInput=[1 1 2 1 1 0]; %[Flag, LQTFP, Contrib, Loutduok,Louthasbo,louthasfs]
        
         axiMeshprep(r,z,n,nang,npanelt,zCoG,depth,w,dir,QTFInput,dirname);% Call the function axiMesh.m
         
          case 9
            %% SOFTWIND PLATFORM
            %% specify n angular 50, npanel 500, z gravity -1.789*Frscaled
        clc
        clear all
        close all
        
        dirname='SOFTWINDNang50Npt1200';
        nang=50;
        npanelt=1200;
        nbfreq=60; % number of calculations = (number of BVP per frequency*number of frequencies)
        
        FrScaled=40;
        R1 =0.14.*FrScaled;
        R2 =0.225.*FrScaled; % Radius of the cylinder
        R121=0.1683*FrScaled;
        R122=0.1966*FrScaled;
        
        Draft =-2.285.*FrScaled; % Height of the submerged part
        
        r=[R1 R1 R121 R122 R2 R2 3*R2/4 2*R2/4 1*R2/4 0]; % r is the first coordinates of the
        z=[0 -0.135.*FrScaled -0.2016*FrScaled -0.2682*FrScaled -0.335.*FrScaled Draft Draft Draft Draft Draft];
        n=length(r);
        zCoG=-1.789*FrScaled;
       
        w= linspace(0.05,3,nbfreq)'; % Periods from 3s to 20s for waves for instance
        dir=0;% angle of the incident waves
        depth=5*FrScaled; % water depth (m)
        QTFInput=[1 1 2 1 1 0]; %[Flag, LQTFP, Contrib, Loutduok,Louthasbo,louthasfs]
        
         axiMeshprep(r,z,n,nang,npanelt,zCoG,depth,w,dir,QTFInput,dirname);% Call the function axiMesh.m
         
          case 10
            %% SOFTWIND PLATFORM
            %% specify n angular 50, npanel 500, z gravity -1.789*Frscaled
        clc
        clear all
        close all
        
        dirname='SOFTWINDNang60Npt1500';
        nang=60;
        npanelt=1500;
        nbfreq=60; % number of calculations = (number of BVP per frequency*number of frequencies)
        
        FrScaled=40;
        R1 =0.14.*FrScaled;
        R2 =0.225.*FrScaled; % Radius of the cylinder
        R121=0.1683*FrScaled;
        R122=0.1966*FrScaled;
        
        Draft =-2.285.*FrScaled; % Height of the submerged part
        
        r=[R1 R1 R121 R122 R2 R2 3*R2/4 2*R2/4 1*R2/4 0]; % r is the first coordinates of the
        z=[0 -0.135.*FrScaled -0.2016*FrScaled -0.2682*FrScaled -0.335.*FrScaled Draft Draft Draft Draft Draft];
        n=length(r);
        zCoG=-1.789*FrScaled;
       
        w= linspace(0.05,3,nbfreq)'; % Periods from 3s to 20s for waves for instance
        dir=0;% angle of the incident waves
        depth=5*FrScaled; % water depth (m)
        QTFInput=[1 1 2 1 1 0]; %[Flag, LQTFP, Contrib, Loutduok,Louthasbo,louthasfs]
        
         axiMeshprep(r,z,n,nang,npanelt,zCoG,depth,w,dir,QTFInput,dirname);% Call the function axiMesh.m
         
           case 11
            %% SOFTWIND PLATFORM
            %% specify n angular 50, npanel 500, z gravity -1.789*Frscaled
        clc
        clear all
        close all
        
        dirname='SOFTWINDNang60Npt2000';
        nang=60;
        npanelt=2000;
        nbfreq=60; % number of calculations = (number of BVP per frequency*number of frequencies)
        
        FrScaled=40;
        R1 =0.14.*FrScaled;
        R2 =0.225.*FrScaled; % Radius of the cylinder
        R121=0.1683*FrScaled;
        R122=0.1966*FrScaled;
        
        Draft =-2.285.*FrScaled; % Height of the submerged part
        
        r=[R1 R1 R121 R122 R2 R2 3*R2/4 2*R2/4 1*R2/4 0]; % r is the first coordinates of the
        z=[0 -0.135.*FrScaled -0.2016*FrScaled -0.2682*FrScaled -0.335.*FrScaled Draft Draft Draft Draft Draft];
        n=length(r);
        zCoG=-1.789*FrScaled;
       
        w= linspace(0.05,3,nbfreq)'; % Periods from 3s to 20s for waves for instance
        dir=0;% angle of the incident waves
        depth=5*FrScaled; % water depth (m)
        QTFInput=[1 1 2 1 1 0]; %[Flag, LQTFP, Contrib, Loutduok,Louthasbo,louthasfs]
        
         axiMeshprep(r,z,n,nang,npanelt,zCoG,depth,w,dir,QTFInput,dirname);% Call the function axiMesh.m
         
             case 12
            %% SOFTWIND PLATFORM
            %% specify n angular 50, npanel 500, z gravity -1.789*Frscaled
        clc
        clear all
        close all
        
        dirname='SOFTWINDNang70Npt2500';
        nang=70;
        npanelt=2500;
        nbfreq=60; % number of calculations = (number of BVP per frequency*number of frequencies)
        
        FrScaled=40;
        R1 =0.14.*FrScaled;
        R2 =0.225.*FrScaled; % Radius of the cylinder
        R121=0.1683*FrScaled;
        R122=0.1966*FrScaled;
        
        Draft =-2.285.*FrScaled; % Height of the submerged part
        
        r=[R1 R1 R121 R122 R2 R2 3*R2/4 2*R2/4 1*R2/4 0]; % r is the first coordinates of the
        z=[0 -0.135.*FrScaled -0.2016*FrScaled -0.2682*FrScaled -0.335.*FrScaled Draft Draft Draft Draft Draft];
        n=length(r);
        zCoG=-1.789*FrScaled;
       
        w= linspace(0.05,3,nbfreq)'; % Periods from 3s to 20s for waves for instance
        dir=0;% angle of the incident waves
        depth=5*FrScaled; % water depth (m)
        QTFInput=[1 1 2 1 1 0]; %[Flag, LQTFP, Contrib, Loutduok,Louthasbo,louthasfs]
        
         axiMeshprep(r,z,n,nang,npanelt,zCoG,depth,w,dir,QTFInput,dirname);% Call the function axiMesh.m
     
end
disp('Please proceed now for the NEMOH preprocessing')


