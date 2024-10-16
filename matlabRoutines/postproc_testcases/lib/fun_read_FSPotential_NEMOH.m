function [w,k,beta,Npanel,Nwline,XX,POT]=fun_read_FSPotential_NEMOH(wd,FileN)
fid=fopen([wd,filesep,'Nemoh.cal'],'r');
for i=1:6
    ligne=fgetl(fid);
end
nBodies=fscanf(fid,'%g',1);
ligne=fgetl(fid);
for i=1:nBodies
    for j=1:3
        ligne=fgetl(fid);
    end
    Ndof=fscanf(fid,'%g',1);
    ligne=fgetl(fid);
    for j=1:Ndof
        ligne=fgetl(fid);
    end
    Nforce=fscanf(fid,'%g',1);
    ligne=fgetl(fid);
    for j=1:Nforce
        ligne=fgetl(fid);
    end
    Nadd=fscanf(fid,'%g',1);
    ligne=fgetl(fid);
    if Nadd>0
        Ndeflm=fscanf(fid,'%g',1);
        ligne=fgetl(fid);
        for j=1:Nadd-1
            ligne=fgetl(fid);
        end
    else
        Ndeflm=0;
    end
end
ligne=fgetl(fid);
ligne=fscanf(fid,'%d',2);
Nw=ligne(2);
ligne=fgetl(fid);
ligne=fscanf(fid,'%g',3);
Nbeta=ligne(1);
betamin=ligne(2);
betamax=ligne(3);
fclose(fid);

fid=fopen([wd,filesep,'QTFPreprocOut',filesep,FileN],'r');
for i=1:4
    ligne=fgetl(fid);
end
ligne=fscanf(fid,'%d',4);
Npanel=ligne(1);
Nwline=ligne(2);
Isym=ligne(3);


ligne=fscanf(fid,'%g',(2^Isym)*3*(Npanel+Nwline));    
XX(1:Npanel,1)=ligne(1:Npanel);
XX(1:Npanel,2)=ligne(Npanel+1:2*Npanel);
XX(1:Npanel,3)=ligne(2*Npanel+1:3*Npanel);
XX(Npanel+1:Npanel+Nwline,1)=ligne(3*Npanel+1:3*Npanel+Nwline);
XX(Npanel+1:Npanel+Nwline,2)=ligne(3*Npanel+Nwline+1:3*Npanel+Nwline*2);
XX(Npanel+1:Npanel+Nwline,3)=ligne(3*Npanel+Nwline*2+1:3*Npanel+Nwline*3);
if Isym==1  
XX(Npanel+Nwline+1:2*(Npanel+Nwline),[1,3])=XX(1:Npanel+Nwline,[1,3]);
XX(Npanel+Nwline+1:2*(Npanel+Nwline),2)=-XX(1:Npanel+Nwline,2);
end
ligne=fgetl(fid);
NPFlow=Npanel+Nwline;
POT=zeros((2^Isym)*NPFlow,2,Nw,Nbeta);
w=zeros(Nw,1);k=w;beta=zeros(Nbeta,1);
for ibeta=1:Nbeta
 for iw=1:Nw
 ligne=fscanf(fid,'%g',3+(2^Isym)*2*(NPFlow));
  w(iw)=ligne(1);
  k(iw)=ligne(2);
  beta(ibeta)=ligne(3);
  POT(1:NPFlow,1,iw,ibeta)=ligne(4:3+NPFlow);
  POT(1:NPFlow,2,iw,ibeta)=ligne(4+NPFlow:3+2*NPFlow);
  if Isym==1  
  POT(NPFlow+1:2*NPFlow,1,iw,ibeta)=ligne(4+2*NPFlow:3+3*NPFlow);
  POT(NPFlow+1:2*NPFlow,2,iw,ibeta)=ligne(4+3*NPFlow:3+4*NPFlow);   
  end
  
  ligne=fgetl(fid);
 end
end



end