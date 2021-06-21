clc
clear all
close all

testcase=10; 

switch testcase 
       case 1
        %  Don't Forget z(i)> z(i+1) !!
        clc
        clear all
        close all
        
        dirname='hemisphereR10InfD';
        nang=50;
        npanelt=600;
        
        
        n1=10; % number of points for the description of the
        n=n1+1; %
        Rayon=10;
        for i=1:n+1
            r(i)=Rayon*cos((n1-i+1)*pi/2/n1-pi/2);
            z(i)=Rayon*sin((n1-i+1)*pi/2/n1-pi/2);
        end
            
        %------------Define calculation options-------------
        nbfreq=20; % number of calculations = (number of BVP per frequency*number of frequencies)
        g=9.811;
        w= linspace(sqrt(g)/nbfreq,sqrt(g),nbfreq)'; % Periods from 3s to 20s for waves for instance
        dir=0;% angle of the incident waves
        depth=600; % water depth (m)
        zCoG=0;
        QTFInput=[0]; %[Flag, LQTFP, Contrib, Loutduok,Louthasbo,louthasfs]
        
        axiMeshprep(r,z,n,nang,npanelt,zCoG,depth,w,dir,QTFInput,dirname);% Call the function axiMesh.m
    case 2
        %  Don't Forget z(i)> z(i+1) !!
        clc
        clear all
        close all
        
        dirname='hemisphereR10InfDLid';
        nang=50;
        npanelt=600;
        
        
        n1=10; % number of points for the description of the
        n=n1+1;
        Rayon=10;
        for i=1:n+1
            r(i)=Rayon*cos((n1-i+1)*pi/2/n1-pi/2);
            z(i)=Rayon*sin((n1-i+1)*pi/2/n1-pi/2);
        end
            rtemp=zeros(1,length(r)+3);ztemp=rtemp;
            rtemp(5:end)=r(1:end-1);ztemp(5:end)=z(1:end-1);
            rtemp(1)=0; rtemp(2)=Rayon/4;rtemp(3)=Rayon/2;rtemp(4)=3*Rayon/4;
            ztemp(1:4)=0;
            r=rtemp;z=ztemp;n=length(r);
        %------------Define calculation options-------------
        nbfreq=20; % number of calculations = (number of BVP per frequency*number of frequencies)
        g=9.811;
        w= linspace(sqrt(g)/nbfreq,sqrt(g),nbfreq)'; % Periods from 3s to 20s for waves for instance
        dir=0;% angle of the incident waves
        depth=600; % water depth (m)
        zCoG=0;
        QTFInput=[1 1 2 1 1 0]; %[Flag, LQTFP, Contrib, Loutduok,Louthasbo,louthasfs]
        
        axiMeshprep(r,z,n,nang,npanelt,zCoG,depth,w,dir,QTFInput,dirname);% Call the function axiMesh.m
      

        case 3
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
         
          case 4
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
         
          case 5
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
         
           case 6
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
         
             case 7
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
         
          case 8
            %% SOFTWIND PLATFORM
            %% specify n angular 50, npanel 500, z gravity -1.789*Frscaled
        clc
        clear all
        close all
        
        dirname='SOFTWINDNang50Npt1000Lid';
        FrScaled=40;
        R1 =0.14.*FrScaled;
        R2 =0.225.*FrScaled; % Radius of the cylinder
        R121=0.1683*FrScaled;
        R122=0.1966*FrScaled;
        
        Draft =-2.285.*FrScaled; % Height of the submerged part
        
        r=[0 R1/4 2*R1/4 3*R1/4 R1 R1 R121 R122 R2 R2 3*R2/4 2*R2/4 1*R2/4 0]; % r is the first coordinates of the
        z=[0 0    0        0  0 -0.135.*FrScaled -0.2016*FrScaled -0.2682*FrScaled -0.335.*FrScaled Draft Draft Draft Draft Draft];
        n=length(r);
        zCoG=-1.789*FrScaled;
        nang=50;
        npanelt=1000;
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
        
        dirname='SOFTWINDNang50Npt1600Lid';
        FrScaled=40;
        R1 =0.14.*FrScaled;
        R2 =0.225.*FrScaled; % Radius of the cylinder
        R121=0.1683*FrScaled;
        R122=0.1966*FrScaled;
        
        Draft =-2.285.*FrScaled; % Height of the submerged part
        
        r=[0 R1/4 2*R1/4 3*R1/4 R1 R1 R121 R122 R2 R2 3*R2/4 2*R2/4 1*R2/4 0]; % r is the first coordinates of the
        z=[0 0    0        0  0 -0.135.*FrScaled -0.2016*FrScaled -0.2682*FrScaled -0.335.*FrScaled Draft Draft Draft Draft Draft];
        n=length(r);
        zCoG=-1.789*FrScaled;
        nang=50;
        npanelt=1600;
        nbfreq=60; % number of calculations = (number of BVP per frequency*number of frequencies)
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
        
        dirname='SOFTWINDNang50Npt2000Lid';
        FrScaled=40;
        R1 =0.14.*FrScaled;
        R2 =0.225.*FrScaled; % Radius of the cylinder
        R121=0.1683*FrScaled;
        R122=0.1966*FrScaled;
        
        Draft =-2.285.*FrScaled; % Height of the submerged part
        
        r=[0 R1/4 2*R1/4 3*R1/4 R1 R1 R121 R122 R2 R2 3*R2/4 2*R2/4 1*R2/4 0]; % r is the first coordinates of the
        z=[0 0    0        0  0 -0.135.*FrScaled -0.2016*FrScaled -0.2682*FrScaled -0.335.*FrScaled Draft Draft Draft Draft Draft];
        n=length(r);
        zCoG=-1.789*FrScaled;
        nang=50;
        npanelt=2000;
        nbfreq=60; % number of calculations = (number of BVP per frequency*number of frequencies)
        w= linspace(0.05,3,nbfreq)'; % Periods from 3s to 20s for waves for instance
        dir=0;% angle of the incident waves
        depth=5*FrScaled; % water depth (m)
        QTFInput=[1 1 2 1 1 0]; %[Flag, LQTFP, Contrib, Loutduok,Louthasbo,louthasfs]
        
         axiMeshprep(r,z,n,nang,npanelt,zCoG,depth,w,dir,QTFInput,dirname);% Call the function axiMesh.m
     
end
disp('Please proceed now for the NEMOH preprocessing')


