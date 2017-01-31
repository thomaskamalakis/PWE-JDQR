function S0=dielectric_averaging2(S);

t0=clock;

% Estimation of the averages of the dielectric constant of the structure
% used in the PWEM calculations

% Inner grid for checking for dielectric interfaces
S.Dxsmall=S.Dx/(S.points_ch-1);                         
S.Dysmall=S.Dy/(S.points_ch-1);
S.Dzsmall=S.Dz/(S.points_ch-1);

% Will store the inverse dielectric constant corresponding to the
% inner grid
S.eintinv=zeros(S.points_ch,S.points_ch,S.points_ch);

% flag used for checking for interfaces
S.mark=zeros(S.Nx,S.Ny,S.Nz);
f=S.iefun;
fin=S.iefunin;
g=@(x,y,z) 1./f(x,y,z);

% Estimate dielectric constant (no averaging)
S.einv=f(S.x,S.y,S.z);

% load points for spherical integration
load -ascii lebedev_011.txt

% assign values to integration angles
theta=lebedev_011(:,1);
phi=lebedev_011(:,2);

% weights for integration
w=lebedev_011(:,3);

% convert angles to radians
theta=theta/180*pi;
phi=phi/180*pi;

% estimate points to be used for the integration in cartesian coordinates
r=min([S.Dx S.Dy S.Dz]);
xi=r*cos(theta).*sin(phi);
yi=r*sin(theta).*sin(phi);
zi=r*cos(phi);

% rsphere=Dx;

vecx=zeros(S.Nx,S.Ny,S.Nz);
vecy=zeros(S.Nx,S.Ny,S.Nz);
vecz=zeros(S.Nx,S.Ny,S.Nz);
S.normvecx=zeros(S.Nx,S.Ny,S.Nz);
S.normvecy=zeros(S.Nx,S.Ny,S.Nz);
S.normvecz=zeros(S.Nx,S.Ny,S.Nz);

S.einv1=zeros(S.Nx,S.Ny,S.Nz);
S.einv2=zeros(S.Nx,S.Ny,S.Nz);
S.e=zeros(S.Nx,S.Ny,S.Nz);


% create coarse voxel to check if there is an interface on nx,ny,nz
xint0=(1:S.points_ch)*S.Dxsmall;
xint0=xint0-mean(xint0);

yint0=(1:S.points_ch)*S.Dysmall;
yint0=yint0-mean(yint0);

zint0=(1:S.points_ch)*S.Dzsmall;
zint0=zint0-mean(zint0);
            
avepoints=9;
Dxave=S.Dx/(avepoints-1);
Dyave=S.Dy/(avepoints-1);
Dzave=S.Dz/(avepoints-1);

xave0=(1:avepoints)*Dxave;
xave0=xave0-mean(xave0);

yave0=(1:avepoints)*Dyave;
yave0=yave0-mean(yave0);

zave0=(1:avepoints)*Dzave;
zave0=zave0-mean(zave0);
              

S.e=g(S.x,S.y,S.z);

% search for interior points

ee=zeros(length(S.x),length(S.y),length(S.z));
for mx=1:length(xint0)
    for my=1:length(yint0)
        for mz=1:length(zint0)            
            ee=ee+f(S.x+xint0(mx),S.y+yint0(my),S.z+zint0(mz));
        end
    end
end
ee=ee/length(xint0)/length(yint0)/length(zint0);

for nx=1:S.Nx
    for ny=1:S.Ny
        for nz=1:S.Nz
            
                         
            % intrface pact
            if ee(nx,ny,nz)~=S.e(nx,ny,nz)
              S.mark(nx,ny,nz)=1;
              % sphere creation to determine vectors
              xx=xi+S.x(nx);
              yy=yi+S.y(ny);
              zz=zi+S.z(nz);
             
              % calculate e for every point on the sphere
              eintsphere=1./fin(xx,yy,zz);
                                  
              if mean(eintsphere)~=min(eintsphere)
                  % perpendicular vector ( where xi*r*eintsphere=r*e)
                  vecx(nx,ny,nz)=4*pi*r^2*sum(w.*xi.*eintsphere);
                  vecy(nx,ny,nz)=4*pi*r^2*sum(w.*yi.*eintsphere);
                  vecz(nx,ny,nz)=4*pi*r^2*sum(w.*zi.*eintsphere);
                  
                  % normalization
                  % set minus normvec. by definition, vec points inside the rod, due to bigger dielectric constant.
                  S.normvecx(nx,ny,nz)= - vecx(nx,ny,nz)/sqrt(vecx(nx,ny,nz).^2+vecy(nx,ny,nz).^2+vecz(nx,ny,nz).^2);
                  S.normvecy(nx,ny,nz)= - vecy(nx,ny,nz)/sqrt(vecx(nx,ny,nz).^2+vecy(nx,ny,nz).^2+vecz(nx,ny,nz).^2);
                  S.normvecz(nx,ny,nz)= - vecz(nx,ny,nz)/sqrt(vecx(nx,ny,nz).^2+vecy(nx,ny,nz).^2+vecz(nx,ny,nz).^2);
              else
                  S.normvecx(nx,ny,nz)= 0.0;
                  S.normvecy(nx,ny,nz)= 0.0;
                  S.normvecz(nx,ny,nz)= 0.0;
              end
              
             % recalculate dielectric in a 6 by 6 voxel, in order to take the average
              
              % create voxel to check if there is an interface on nx,ny,nz
              xave=xave0+S.x(nx);
              yave=yave0+S.y(ny);
              zave=zave0+S.z(nz);

              % calculate e for every point in the voxel
              eintave_inv=f(xave,yave,zave);
              eintave=1./eintave_inv;
              
              S.einv1(nx,ny,nz)=mean(mean(mean(1./(eintave))));
              S.einv2(nx,ny,nz)=1/mean(mean(mean(eintave)));
            else
                S.einv1(nx,ny,nz)=S.einv(nx,ny,nz);
                S.einv2(nx,ny,nz)=S.einv(nx,ny,nz);
            end
        end
    end
    
end
t_cur=clock;
elapsed_time=etime(t_cur,t0);
S.avgtime=elapsed_time;
S.einv0=mean(mean(mean(S.einv)));
S0=S;
