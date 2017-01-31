
function [vec_a,vec_b, kplus_M]=orth_vecs(k,Gx,Gy,Gz);

Nx=length(Gx);
Ny=length(Gy);
Nz=length(Gz);

vec_a=zeros(Nx,Ny,Nz,3);
vec_b=zeros(Nx,Ny,Nz,3);
kplus_M=zeros(Nx,Ny,Nz);

for cx=1:Nx
    for cy=1:Ny
        for cz=1:Nz

            Gc=[Gx(cx) Gy(cy) Gz(cz)];
            kplusG=[k(1)+Gc(1) k(2)+Gc(2) k(3)+Gc(3)];
            kplusGt=[k(1)+Gc(1) k(2)+Gc(2) 0];
            
            kplus_M(cx,cy,cz)=norm(kplusG);
            
            if ( Gc(3) > 0 )
               tmpa=[ -(k(2)+Gc(2)) (k(1)+Gc(1)) 0 ]/norm(kplusGt);
               tmpb=1/norm(kplusG)*[ -(k(1)+Gc(1))*Gc(3)/norm(kplusGt) -(k(2)+Gc(2))*Gc(3)/norm(kplusGt) norm(kplusGt)];

            elseif ( Gc(3) == 0 )
               tmpa=1/sqrt(2)/norm(kplusGt)*[k(2)+Gc(2) -(k(1)+Gc(1)) norm(kplusGt)];
               tmpb=1/sqrt(2)/norm(kplusGt)*[k(2)+Gc(2) -(k(1)+Gc(1)) -norm(kplusGt)];

            elseif ( Gc(3) < 0 )
               tmpb=[ -(k(2)+Gc(2)) (k(1)+Gc(1)) 0 ]/norm(kplusGt);
               tmpa=1/norm(kplusG)*[ (k(1)+Gc(1))*Gc(3)/norm(kplusGt) (k(2)+Gc(2))*Gc(3)/norm(kplusGt) -norm(kplusGt)];
            end

            vec_a(cx,cy,cz,:)=tmpa;
            vec_b(cx,cy,cz,:)=tmpb;
        end 
    end
end

