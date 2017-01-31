%
% COMPARISON OF INTERIOR POINT EIGENVALUE ESTIMATION METHODS
%

% Clear everything
clc
clear all
close all
tic
% Global variable measuring the number of calls
global nocalls
nocalls=0;

% Constants
S.e0=8.854e-12;
S.m0=4*pi*1e-7;

% supercell size
S.a=1;          % Cell size along x (in normalized units)
S.b=10.5*S.a;   % Cell size along y (in normalized units)
S.c=4*S.a;      % Cell size along z (in normalized units)
S.a0=440e-9;    % lattice constant

% Define supercell geometry
S.ra=0.2862*S.a;           % hole radii (normalized by the lattice constant)
S.eup=1;                   % dielectric constant above the slab
S.edown=1;                 % dielectric constant above the slab
S.eslab=12.00;             % dielectric constant of slab
S.ehole=1;                 % dielectric constant of the holes
S.h=0.5314*S.a;            % slab thickness

% rod positions
S.rodsx=S.a*[0 1/2 -1/2 0 1/2 -1/2 0 1/2 -1/2 ...
         1/2 -1/2 0 1/2 -1/2 0 1/2 -1/2 0];

S.rodsy=S.a*[-3*sqrt(3), -5*sqrt(3)/2, -5*sqrt(3)/2, -2*sqrt(3), ...
         -3*sqrt(3)/2, -3*sqrt(3)/2, -sqrt(3), -sqrt(3)/2, -sqrt(3)/2 ...
         sqrt(3)/2, sqrt(3)/2, sqrt(3) 3*sqrt(3)/2 3*sqrt(3)/2 ...
         2*sqrt(3) 5*sqrt(3)/2 5*sqrt(3)/2 3*sqrt(3)];

S.rodsr=S.ra*ones(size(S.rodsy));

% number of plane waves
Nplws=8;
  
S.Nplw=Nplws;
S.Nx=S.Nplw;            % Number of plane waves along x
S.Ny=4*S.Nplw;          % Number of plane waves along y
S.Nz=2*S.Nplw;          % Number of plane waves along z
S.N=S.Nx*S.Ny*S.Nz;     % Total number of plane waves

% Number of bands
S.Nb=1;

% Averaging parameters
S.points_ch=3;                    % use a 3x3 cube to determine whether there is a boundary at one single point

% Build lattice and reciprocal vectors 
S=lattice_vectors_rect(S);

% Set the dielectric constant distirbution to be that of an air membrane
S.iefun=@(x,y,z) slab_holes(S.eslab,S.ehole,S.eup,S.edown,x,y,z,S.rodsx,S.rodsy,S.rodsr,S.h);
S.iefunin=@(x,y,z) slab_holes2(S.eslab,S.ehole,S.eup,S.edown,x,y,z,S.rodsx,S.rodsy,S.rodsr,S.h);

fprintf(1,'Starting calculations for %d plane waves\n',S.Nplw);

% perform dielectric averaging
S=dielectric_averaging2(S);

% calculate dielectric tensors
S=dielectric_tensors(S);
fprintf(1,'Completed averaging in %6.1f secs\n',S.avgtime);

% area to calculate confinement factor
S.czmin=-S.h/2;                              % bottom of slab
S.czmax=+S.h/2;                              % top of slab
S.cymin=-S.a;                                % one rod away to the left
S.cymax=+S.a;                                % one rod away to the left
S.cxmin=-S.a/2;                              % start of supercell along propagation direction
S.cxmax=+S.a/2;                              % end of supercell along propagation direction

% wavevector in question
S.k=[1.0*pi/S.a 0 0];

% Calculation of the orthogonal vectors
S=st_orth_vecs(S.k,S);

% Band-type
S.opt.bandtype='all';

% Calculate 
S=set_operator2(S);

wca_t=0.28;
tau=(2*pi*S.a*wca_t)^2;
S.solvopts.tau=tau;
S.solvopts.tol=1e-6;
S.solvopts.v0=randn(2*S.Nx*S.Ny*S.Nz,1);
S.solvopts.n=2*S.Nx*S.Ny*S.Nz;

% search space size after implicit restart   
if ~isfield(S.solvopts,'jmin')
   S.solvopts.jmin=20;
end

% search space size before implicit restart
if ~isfield(S.solvopts,'jmax')
   S.solvopts.jmax=S.solvopts.jmin+20;
end

% target eigenvalue
if ~isfield(S.solvopts,'tau')
   S.solvopts.tau=0.0;
end

% number of eigenvalues
if ~isfield(S.solvopts,'neig')
   S.solvopts.neig=1;
end

% initial start point
if ~isfield(S.solvopts,'v0')
   S.solvopts.v0=randn(S.solvopts.n,1);
end

% residual tolerance
if ~isfield(S.solvopts,'tol')
  S.solvopts.tol=1e-6;
end 

% solver type for the correction equation.
if ~isfield(S.solvopts,'solver')
   S.solvopts.solver='SIMPLE';
end

% max number of iterations       
if ~isfield(S.solvopts,'maxiters')
   S.solvopts.maxiters=100;
end

% for GMRES check additional fields.
if strcmp(S.solvopts.solver,'GMRES')

   if ~isfield(S.solvopts,'rstgmres')
       S.solvopts.rstgmres=5;
   end

   if ~isfield(S.solvopts,'maxitersgmres')
       S.solvopts.maxitersgmres=S.solvopts.rstgmres;
   end

   if ~isfield(S.solvopts,'tolgmres')
       S.solvopts.tolgmres=1e-6;
   end

end

% Preconditionner    
if ~isfield(S.solvopts,'K')
   S.solvopts.K=@(V,l) PVappx(V,S,l);
end

% Verbose level
if ~isfield(S.solvopts,'vrblevel')
   S.solvopts.vrblevel=1;
end

% operator
S.solvopts.funA=@(V) MV(V,S);
    
options=S.solvopts;
funA=S.solvopts.funA;
[v,lambda,output]=jdqr(funA,options);
S.solvoutput=output;
S.v=v;
S.lambda=lambda;
toc
