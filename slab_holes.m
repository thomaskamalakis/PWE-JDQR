function einv=slab_holes(eslab,ehole,eup,edown,x,y,z,rodsx,rodsy,rodsr,h)

% einv=inv_dielectric_pcslab(ea,eb,ec,x,y,z,rodsx,rodsy,rodsr,h);
%
% returns the inverse dielectric constant of a photonic crystal slab of
% thickness h.
% The coordinates of the centers of the rods in the (x,y) plane are
% described by rodsx,rodsy while rodsr contains the radii of the rods
% eup is the dielectric constant of the rods above the slab
% edown is the dielectric constant of the rods below the slab
% eslab is the dielectric constant of the slab
% ehole is the dielectric constant of the holes
% x, y and z are the grid axis

Nx=length(x);
Ny=length(y);
Nz=length(z);
Nrods=length(rodsx);

x=reshape(x,Nx,1,1);
y=reshape(y,1,Ny,1);
z=reshape(z,1,1,Nz);

xxx=x(:,ones(1,Ny),ones(1,Nz));
yyy=y(ones(1,Nx),:,ones(1,Nz));
zzz=z(ones(1,Nx),ones(1,Ny),:);

up=abs(zzz>=h/2);
down=abs(zzz<=-h/2);
inside=zeros(size(zzz));
for m=1:Nrods
    inside0=((xxx-rodsx(m)).^2 + (yyy-rodsy(m)).^2 <= rodsr(m)^2) & ( abs(zzz) < h/2 );
    inside=inside | inside0;
end
slab=~up & ~down & ~inside;
inside=abs(inside);
slab=abs(slab);

einv = 1/ehole * inside + 1/eslab * slab + 1/eup * up + 1/edown * down;











