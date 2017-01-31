function Vout=MV(V,S)

global nocalls DMVtime
t1=clock;

if exist('nocalls')==0
    nocalls=0;
else
    nocalls=nocalls+1;
end
   
fftshift3=@(x) fftshift(fftshift(fftshift(x,1),2),3);

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

% Now you need to reshape in order to get 3D matrices and work on the
% FFT/IFFT

% Apply the first curl
Hp2=1j*Hq1.*S.kplus_M;
Hq2=-1j*Hp1.*S.kplus_M;

% Project on the x,y,z axes
Hx=Hp2.*S.vec_a(:,:,:,1)+Hq2.*S.vec_b(:,:,:,1);
Hy=Hp2.*S.vec_a(:,:,:,2)+Hq2.*S.vec_b(:,:,:,2);
Hz=Hp2.*S.vec_a(:,:,:,3)+Hq2.*S.vec_b(:,:,:,3);

% Calculate field in the space domain
Hx2=fftshift3(fftn(fftshift3(Hx)));
Hy2=fftshift3(fftn(fftshift3(Hy)));
Hz2=fftshift3(fftn(fftshift3(Hz)));

% Multiply with the inverse of the dielectric constant
Hx3=S.Txx.*Hx2+S.Txy.*Hy2+S.Txz.*Hz2;
Hy3=S.Txy.*Hx2+S.Tyy.*Hy2+S.Tyz.*Hz2;
Hz3=S.Txz.*Hx2+S.Tyz.*Hy2+S.Tzz.*Hz2;

% Bring the field back in the plane wave basis
Hx3=fftshift3(ifftn(fftshift3(Hx3)));
Hy3=fftshift3(ifftn(fftshift3(Hy3)));
Hz3=fftshift3(ifftn(fftshift3(Hz3)));

% Project along the unitary vectors vec_a and vec_b
Hp3=Hx3.*S.vec_a(:,:,:,1)+Hy3.*S.vec_a(:,:,:,2)+Hz3.*S.vec_a(:,:,:,3);
Hq3=Hx3.*S.vec_b(:,:,:,1)+Hy3.*S.vec_b(:,:,:,2)+Hz3.*S.vec_b(:,:,:,3);

% Calculate the second curl
Hp4=1j*Hq3.*S.kplus_M;
Hq4=-1j*Hp3.*S.kplus_M;

% Finally reshape the plane waves into on big happy vector
Hp5=reshape(Hp4,N,1);
Hq5=reshape(Hq4,N,1);
Vout=[Hp5;Hq5];

t2=clock;
if isempty(DMVtime)
    DMVtime=etime(t2,t1);
else
    DMVtime=DMVtime+etime(t2,t1);
end


