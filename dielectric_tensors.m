function S0=dielectric_tensors(S);

S.Txx=zeros(S.Nx,S.Ny,S.Nz);
S.Tyy=zeros(S.Nx,S.Ny,S.Nz);
S.Tzz=zeros(S.Nx,S.Ny,S.Nz);
S.Txy=zeros(S.Nx,S.Ny,S.Nz);
S.Txz=zeros(S.Nx,S.Ny,S.Nz);
S.Tyz=zeros(S.Nx,S.Ny,S.Nz);

for nx=1:S.Nx
    for ny=1:S.Ny
        for nz=1:S.Nz
            if S.mark(nx,ny,nz)==0
                Txx=S.einv(nx,ny,nz);
                Tyy=S.einv(nx,ny,nz);
                Tzz=S.einv(nx,ny,nz);
                Txy=0;
                Txz=0;
                Tyz=0;
            else
                Txx=S.einv1(nx,ny,nz)*(S.normvecx(nx,ny,nz)*S.normvecx(nx,ny,nz)) + S.einv2(nx,ny,nz)*(1-S.normvecx(nx,ny,nz)*S.normvecx(nx,ny,nz));
                Tyy=S.einv1(nx,ny,nz)*(S.normvecy(nx,ny,nz)*S.normvecy(nx,ny,nz)) + S.einv2(nx,ny,nz)*(1-S.normvecy(nx,ny,nz)*S.normvecy(nx,ny,nz));
                Tzz=S.einv1(nx,ny,nz)*(S.normvecz(nx,ny,nz)*S.normvecz(nx,ny,nz)) + S.einv2(nx,ny,nz)*(1-S.normvecz(nx,ny,nz)*S.normvecz(nx,ny,nz));
                Txy=S.einv1(nx,ny,nz)*(S.normvecx(nx,ny,nz)*S.normvecy(nx,ny,nz)) + S.einv2(nx,ny,nz)*(-S.normvecx(nx,ny,nz)*S.normvecy(nx,ny,nz));
                Txz=S.einv1(nx,ny,nz)*(S.normvecx(nx,ny,nz)*S.normvecz(nx,ny,nz)) + S.einv2(nx,ny,nz)*(-S.normvecx(nx,ny,nz)*S.normvecz(nx,ny,nz));
                Tyz=S.einv1(nx,ny,nz)*(S.normvecy(nx,ny,nz)*S.normvecz(nx,ny,nz)) + S.einv2(nx,ny,nz)*(-S.normvecy(nx,ny,nz)*S.normvecz(nx,ny,nz));
            end
            S.Txx(nx,ny,nz)=Txx;
            S.Txy(nx,ny,nz)=Txy;
            S.Txz(nx,ny,nz)=Txz;
            
            S.Tyx(nx,ny,nz)=Txy;
            S.Tyy(nx,ny,nz)=Tyy;
            S.Tyz(nx,ny,nz)=Tyz;
            
            S.Tzx(nx,ny,nz)=Txz;
            S.Tzy(nx,ny,nz)=Tyz;
            S.Tzz(nx,ny,nz)=Tzz;
                        
            T=[Txx Txy Txz; Txy Tyy Tyz; Txz Tyz Tzz];
            iT=inv(T);
            iTxx=iT(1,1);
            iTxy=iT(1,2);
            iTxz=iT(1,3);
            iTyx=iT(2,1);
            iTyy=iT(2,2);
            iTyz=iT(2,3);
            iTzx=iT(3,1);
            iTzy=iT(3,2);
            iTzz=iT(3,3);
                        
            S.iTxx(nx,ny,nz)=iTxx;
            S.iTxy(nx,ny,nz)=iTxy;
            S.iTxz(nx,ny,nz)=iTxz;
            
            S.iTyx(nx,ny,nz)=iTyx;
            S.iTyy(nx,ny,nz)=iTyy;
            S.iTyz(nx,ny,nz)=iTyz;
            
            S.iTzx(nx,ny,nz)=iTzx;
            S.iTzy(nx,ny,nz)=iTzy;
            S.iTzz(nx,ny,nz)=iTzz;
            
        end
    end
end
t_cur=clock;

S0=S;