function einv=slab_holes2(eslab,ehole,eup,edown,x,y,z,rodsx,rodsy,rodsr,h)

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

Nrods=length(rodsx);

up=abs(z>=h/2);
down=abs(z<=-h/2);
inside=zeros(size(z));
for m=1:Nrods
    inside0=((x-rodsx(m)).^2 + (y-rodsy(m)).^2 <= rodsr(m)^2) & ( abs(z) < h/2 );
    inside=inside | inside0;
end
slab=~up & ~down & ~inside;
inside=abs(inside);
slab=abs(slab);

einv = 1/ehole * inside + 1/eslab * slab + 1/eup * up + 1/edown * down;











