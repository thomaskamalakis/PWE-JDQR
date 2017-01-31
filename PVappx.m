function Vout=PVappx(V,S,tau)
% Diagonal approximation for MV

% First find out the plane waves
Nx=S.Nx;
Ny=S.Ny;
Nz=S.Nz;

% Then take out the p and the q components of the magnetic field plane
% waves
N=Nx*Ny*Nz;
Hp1=V(1:N);
hp1=reshape(Hp1,Nx,Ny,Nz);
Hp1=hp1;
Hq1=V(N+1:2*N);
hq1=reshape(Hq1,Nx,Ny,Nz);
Hq1=hq1;
einv0=S.einv0;

Hp1=1/einv0./(S.kplus_M.^2-einv0.*tau).*Hp1;
Hq1=1/einv0./(S.kplus_M.^2-einv0.*tau).*Hq1;

Hp1=reshape(Hp1,N,1);
Hq1=reshape(Hq1,N,1);

Vout=[Hp1;Hq1];


