function S0=lattice_vectors_rect(S);

% build coordinate axis
mx=-S.Nx/2:1:(S.Nx/2-1);
my=-S.Ny/2:1:(S.Ny/2-1);
mz=-S.Nz/2:1:(S.Nz/2-1);

% Build reciprocal lattice vectors
S.Dx=S.a/S.Nx;
S.Dy=S.b/S.Ny;
S.Dz=S.c/S.Nz;

% points in the x,y,z axis
S.x=mx*S.Dx;
S.y=my*S.Dy;
S.z=mz*S.Dz;

% Build lattice vectors
S.DGx=2*pi/S.a;
S.DGy=2*pi/S.b;
S.DGz=2*pi/S.c;

S.Gx=mx*S.DGx;
S.Gy=my*S.DGy;
S.Gz=mz*S.DGz;

S0=S;
