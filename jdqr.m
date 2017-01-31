function [v,lambda,output]=jdqr(A,options)

% 
% [v,lambda]=jdqr(A,options);
%
% calculation of a few interior eigenvalues of the Hermitian matrix A using 
% Jacobi-Davidson style QR algorithm. The parameters of the methods are
% determined in the structure options.
%
% option.jmin is the dimension of the search space after an implicit
%             restart. If unspecified we assume jmin=5
%
% option.jmax is the dimension of the search space before an implicit
%             restart. If unspecified we assume jmax=jmin+5
%
% option.tau is the target for the eigenvalues. The algorithm searches for
%            interior eigenvalues near this target value. If unspecified we
%            assume tau=0
% option.maxiters is the maximum number of restarts allowed (=100 if left
%            unspecified)
%
% options.n is the size of the matrix
%
% options.neig is the number of eigenvalues that need to be calculated. If
%              unspecified the algorithm computes one eigenvalue
%
% options.v0 is a n x 1 starting vector for the eigenproblem solution where
%            n=options.n
%
% options.tol is the tolerance of the method.
%
% options.K is the preconditionner matrix. If left unspecified, no
%           preconditionning is used.
%
% options.solver determines the type of solver used for the correction
%                equation: options.solver='SIMPLE' implements a simple step
%                approximate solution, options.solver='GMRES' uses the
%                GMRES method to solve the correction equation. Defaults to
%                'SIMPLE'
%
% options.rstgmres is the size of the GMRES space used in solving the
%                  correction equation before restarting GMRES (=5 if
%                  unspecified)
% options.maxitersgmres is the maximum number of GMRES iterations
%                  (=options.rstgmres if unspecified)
% options.tolgmres is the specified tolerance for GMRES termination (1e-6
%                  if left unspecified).
% options.verblvl is the verbose level. If equal to 0 nothing is shown
%                 during the execution. If equal to 1, the consecutive
%                 residues and eigenvalue estimates are displayed
%
% A and options.K may be function handles.
%
% output returns more detailed output of the method
%
% output.Q and output.R are partial Schur decompositions of A calculated by
%                       the algorithm. 
% output.res is the residue for the last iteration
% output.inner is the number of inner iterations performed
% output.outer is the number of outer iterations performed
% output.reshis contains the history of the residues in the inner
%               iterations of the algorithm
% output resindx is simply a vector from 1 to output.inner that can be used
%                to plot the output data
% output.lamhis contains the history of eigenvalue computations in the
%               interior of the algorithm

global nocalls

t0=clock;

if isfield(options,'jmin')
    jmin=options.jmin;
else
    jmin=5;
end

if isfield(options,'jmax')
    jmax=options.jmax;
else
    jmax=jmin+5;
end

if isfield(options,'neig')
    kmax=options.neig;
else
    kmax=1;
end
    
if isfield(options,'n')
    n=options.n;
else
    fprintf('Error: must specify size of matrix in options structure\n');
    return;
end

if isfield(options,'v0')
    v0=options.v0;
else
    v0=randn(n,1);    
end

if isfield(options,'tol')
    tol=options.tol;
else
    tol=1e-6;
end

if isfield(options,'maxiters')
    maxiters=options.maxiters;
else   
    maxiters=100;
end

if isfield(options,'K')
    K=options.K;
else
    K=@(x) x;
end

if isfield(options,'tau')
    tau=options.tau;
else
    tau=0.0;
end

if isfield(options,'solver')
    solver=options.solver;
else
    solver='SIMPLE';
end

if strcmp(solver,'GMRES')
   
    if isfield(options,'Agmres')
        Agmres=options.Agmres;
    else
        Agmres=A;
    end

    if isfield(options,'rstgmres')
        restartgmres=options.rstgmres;
    else
        restartgmres=5;
    end

    if isfield(options,'tolgmres')
        tolgmres=options.tolgmres;
    else
        tolgmres=1e-6;
    end

    if isfield(options,'maxitersgmres')
        maxitersgmres=options.maxitersgmres;
    else
        maxitersgmres=restartgmres;
    end
    
    if isa(Agmres,'double')
        funAgmres=@(x) Agmres*x;
    else
        funAgmres=Agmres;
    end
end

if isfield(options,'vrblevel')
    vrb=options.vrblevel;
else
    vrb=1;
end

% Check whether matrices are function handles or actual matrices

if isa(A,'double')
    funA=@(x) A*x;
else
    funA=A;
end


if isa(K,'double')
    funK=@(x) K*x;
else
    funK=K;
end


% Initializations
Q=[];
R=[];
Y=[];
H=[];
V=[];
VA=[];
M=[];
k=0;
j=0;
nor=0;
inner=0;
noresiduals=jmax+maxiters*(jmax-jmin);
reshis=zeros(1,noresiduals);
lamhis=zeros(1,noresiduals);
timehis=zeros(1,noresiduals);
callshis=zeros(1,noresiduals);

while (k<kmax) && (nor<maxiters)

  inner=inner+1;
  
  if (j==0)
    v=v0;
  else
    
    % Calculate RHS
    if nargin(funK)==1
        rold=r;
        r=funK(r)-Ytilde/Htilde*(Qtilde'*funK(r));
    else
        rold=r;
        r=funK(r,lambda)-Ytilde/Htilde*(Qtilde'*funK(r,lambda));
    end
       
    % approximate solution of the correction equation
    
    if strcmp(solver,'SIMPLE')
       v=-funK(rold,lambda);       
       v=projQ(Q,v);
    elseif strcmp(solver,'GMRES')
        invHtilde=inv(Htilde);
        [v,itersgmres,res]=gmresjdqr(funAgmres,r,Ytilde,invHtilde,Qtilde,funK,lambda,1e-9*ones(n,1),restartgmres,tolgmres,maxitersgmres);
    end
    
   end
  
  % projected problem
  
  v=rmgs(V,v);
  v=v/norm(v);
  vA=funA(v);
  
  if isempty(V)
    M=v'*vA;
  else
    M=[M, V'*vA; v'*VA, v'*vA];
  end
    
  V=[V v];
  VA=[VA vA];
  [U,S]=schur(M);
  
  % Schur form sorting
  [U,S]=schursort(U,S,tau);
  
  j=j+1;
  found=true;
  
  while found
      
      % Ritz approximation
      lambda=S(1,1);     
      q=V*U(:,1);
      if nargin(funK)==1         
        y=funK(q);
      else
        y=funK(q,lambda);
      end
      r=funA(q)-lambda*q;
      [r,s,~]=rmgs(Q,r);
      Qtilde=[Q q];
      Ytilde=[Y y];
          
      reshis(inner)=norm(r);
      lamhis(inner)=lambda;
      timehis(inner)=etime(clock,t0);
      callshis(inner)=nocalls;
      if vrb
          if exist('nocalls')==0
              fprintf(1,'Inner iteration %d, eigenvalue estimate:%6.2f, residue %e\n',inner,lambda,reshis(inner));
          else               
              fprintf(1,'Inner iteration %d, eigenvalue estimate:%6.2f, residue %e, no calls=%d\n',inner,lambda,reshis(inner),nocalls);
          end
          
      end
      
      if isempty(H)
        Htilde=q'*y;
      else
        Htilde=[H, Q'*y; q'*Y q'*y];
      end
            
      % found and implicit restart
      
      found=( norm(r)<=tol ) & (j>1 | k==kmax-1);
      
      if found
          
          Q=Qtilde;
          R=[R, s; zeros(1,k), lambda];
          k=k+1;
          
          if k==kmax
              break;
          end
          
          Y=Ytilde;
          H=Htilde;
          J=2:j;
          j=j-1;
          
          V=V*U(:,J);
          VA=VA*U(:,J);
          S=S(J,J);
          M=S;
          U=eye(size(M));
          
      elseif (j==jmax)
      
          if vrb
            fprintf(1,'\nOuter iteration %d\n',nor);
          end
                    
          nor=nor+1;          
          j=jmin;
          J=[1:j];
          V=V*U(:,J);
          VA=VA*U(:,J);
          S=S(J,J);
          M=S;
          U=eye(size(M));
      end
      
  end
  
end
output.R=R;
output.Q=Q;
output.inner=inner;
output.outer=nor;
output.res=norm(r);
output.reshis=reshis(1:inner);
output.resindx=1:inner;
output.lamhis=lamhis(1:inner);
output.timehis=timehis(1:inner);
output.callshis=callshis(1:inner);
R=R(2:end,:);
lambda=diag(R);
v=Q;

end

function [x,iters,r]=gmresjdqr(funA,r,Ytilde,invHtilde,Qtilde,funK,lambda,x0,m,tol,maxiters)

j=1;
Y=Ytilde;
Q=Qtilde;
Hi=invHtilde;
Aproj1=@(x) x-Y*Hi*(Q'*x);
if nargin(funK)==2
    funK2=funK;
    funK=@(x) funK2(x,lambda);
end
Aproj2=@(x) funK(funA(x)-lambda*x);
Ad=@(x) Aproj1(Aproj2(Aproj1(x)));

% right hand side
f=-r;

r=norm(f-Ad(x0))/norm(f);
xm=x0;

while (r>=tol) && (j<=maxiters)   

    r0=f-Ad(x0);
    beta=norm(r0);
    v1=r0/beta;

    % Arnoldi factorization
    [V,~,~,~,Haug]=arnoldid(Ad,m,v1);

    % QR factorization
    [Q1,R]=qr(Haug);
    Q1=Q1';

    % Calculate minimization problem
    e1=zeros(m+1,1);
    e1(1)=1;
    gm=beta*Q1*e1;

    % remove last line of R
    R=R(1:m,:);

    % remove last element of gk
    gm=gm(1:m);

    % calculate solution of minimization problem
    ym=R\gm;    
    xm=x0+V*ym;
    
    % project to a space orthogonal to Q
    xm=projQ(Q,xm);
    
    x0=xm;
    r=norm(f-Ad(xm))/norm(f);    
    %fprintf(1,'GMRES iteration %d, residue=%e\n',j,r);
    j=j+1;
end

x=xm;
iters=j-1;
end

function [V,H,hk_1,vk_1,Haug]=arnoldid(funA,k,v)
%
% [V,H,hk_1,vk_1,Haug]=arnoldid(funA,k,v)
%
% Arnoldi factorization for the square matrix A with size(A)=[n,n].
% Returns a Hessenber matrix H with size(H)=[k,k] and a unitary matrix V with size(V)=[n,k], such that
% V'*V=eye(k,k). The matrices are selected so that
% A*V=V*H+hk_1*v_k1*ek', where size(ek)=[k,1], ek(k)=1 and zero otherwise
% Note that V'*A*V=H
% Haug is the augmented H matrix, Haug=[H; [0 ... 0 hk_1];
% funA(x) is a function returning A*x
% k is the order of the factorization
% v is the initial vector
% 

ck=1;
H=zeros(k,k);
v1=v/norm(v);
n=size(v,1);
V=zeros(n,k);
vknew=v1;

while ck<=k
    V(:,ck)=vknew;
    w=funA(V(:,ck));
    for i=1:ck        
        H(i,ck)=V(:,i)'*w;
    end   
    for i=1:ck
        w=w-H(i,ck)*V(:,i);
    end
    hknew=norm(w);
    if ck+1<=k
        H(ck+1,ck)=hknew;
    end    
    vknew=w/hknew;      
    ck=ck+1;
end
hk_1=hknew;
vk_1=vknew;

Haug=[H; [zeros(1,k-1) hk_1]];
    
end

function [U1,D1]=schursort(U,D,tau)
%
% [U1,D1]=schursort(U,D,tau);
%
% Sorts Schur decomposition so that abs(d(i)-tau) >= abs (d(i-1)-tau)
% where d(i)=D(i,i)
%

d=diag(D);
[~,i]=sort(abs(d-tau));
D1=d(i);
D1=diag(D1);
U1=U(:,i);
      
end

function [v,y,vn]=rmgs(V,w)

[m,n]=size(V);
v=w;

% orthogonalize w with respect to V

for k=1:n
   v=v-(V(:,k)'*w)*V(:,k);      
end
vn=v/norm(v);
V=[V vn];
y=V'*w;
end

function z=projQ(Q,y);

Q=mymgs(Q);
[n,m]=size(Q);
z=y;

for k=1:m
    z=z-(Q(:,k)'*z)*Q(:,k);
end

end


    
    

