% Node-based topology optimization 
function [l,v]=sample3

% Problem settings
top.nx = 120;   % The number of elements in the horizontal axis
top.ny = 80;   % The number of elements in the vertical axis
top.vol  = .3;  % Volume constraint
top.pnl  = 3;   % Penalization parameter


% Parameters for optimization
opt.sw1=0;  % OC
%opt.sw1=1; % SLP
opt.convprm1=0.01;
opt.convprm2=20;
% Initialization
rho(1:top.ny+1,1:top.nx+1) = top.vol;
%Edist=zeros(top.ny,top.nx);
itr = 0;
change = 1.;
minc=1.e10;
convflag=0;

% Sensitiity for volume constraint
dvdr=ones(top.ny+1,top.nx+1);
dvdr(1,:)=dvdr(1,:)*0.5;
dvdr(:,1)=dvdr(:,1)*0.5;
dvdr(top.ny+1,:)=dvdr(top.ny+1,:)*0.5;
dvdr(:,top.nx+1)=dvdr(:,top.nx+1)*0.5;

% Optimization loop begins here
while change > 0.01 && convflag <=20
    itr = itr + 1;
    rhoold = rho;
    [U,l,Edist]=fea(top,rho);
    [dldr]=fesens(U,top,rho);

    if opt.sw1==0 % Design variable update
        [rho] = OC(top,rho,dldr,dvdr); % OC
    else
        [rho] = SLP(top,rho,dldr,dvdr); % SLP
    end

    % Convergence check
    change = max(max(abs(rho-rhoold)));
    if minc>l
        minc=l;
        convflag=0;
    else
        convflag=convflag+1;
    end

    v=sum(sum(rho.*dvdr)) / (top.nx*top.ny);
    % Show optimization process data
    fprintf('It.:%4i Cmp.:%8.3f Chng.:%7.3f Vol.:%6.3f\n', itr, l, change, v);
    % Show density distribution
    colormap(gray);imagesc(1-Edist);caxis([0 1]);
    axis equal;axis off;drawnow;

end
end

%%%%% Optimality Creteria %%%%%
function [xnew]=OC(top,rho,dldr,dvdr)
ny=top.ny;
nx=top.nx;
vol=top.vol;
lambda1 = 0; lambda2 = 1e4;
mvlmt = 0.1; %Move-limit
eta=0.5;  %Damping factor

A=-dldr./dvdr;
while (lambda2-lambda1 > 1e-4)
  lmid = 0.5*(lambda2+lambda1);
  xnew = max(0.00,max(rho-mvlmt,min(1.,min(rho+mvlmt,rho.*(A/lmid).^eta))));
  if sum(sum(xnew.*dvdr)) - vol*nx*ny > 0;
    lambda1 = lmid;
  else
    lambda2 = lmid;
  end
end
end

%%%%% Sequential Linear Programming %%%%%
function [rho]=SLP(top,rho,dldr,dvdr)
ny=top.ny;
nx=top.nx;
vol=top.vol;
nvar=(nx+1)*(ny+1); %number of variables
mvlmt = 0.1; %Move-limit
LB = 0.00;    %Lower bound
UB = 1.;       %Upper bound
rho=reshape(rho,nvar,1);

%side constraints
LBm=max(LB,rho*(1-mvlmt));
UBm=min(UB,rho*(1+mvlmt));

%volume constraint
b(1)=vol*nx*ny;
A=reshape(dvdr,1,nvar);

%Solve linear programming problem
rho = reshape(linprog(dldr,A,b,[],[],LBm,UBm),ny+1,nx+1);
end

%%%%% Finite element anlysis %%%%%
function [U, l, Edist]=fea(top,rho)
ny=top.ny;
nx=top.nx;
pnl=top.pnl;
K=zeros(2*(nx+1)*(ny+1), 2*(nx+1)*(ny+1));
F=sparse(2*(ny+1)*(nx+1),1);
U=sparse(2*(ny+1)*(nx+1),1);
Edist=zeros(ny,nx);
for iy = 1:ny
for ix = 1:nx
  n1=(ny+1)*(ix-1)+iy;
  n2=(ny+1)* ix +iy;
  elem=[2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1;2*n2+2;2*n1+1; 2*n1+2];
  enode=[n1; n2; n2+1; n1+1] ;
  xc=[ix-1 ix   ix ix-1];
  yc=[iy-1 iy-1 iy iy];
  [KE, EE]=camd(rho(enode),pnl,xc,yc);
  K(elem,elem)=K(elem,elem)+KE;
  Edist(iy,ix)=EE;
end
end

% Boundary conditions
% Short cantilever
F(2*(nx*(ny+1)+(ny+2)/2-2),1)=-1;
F(2*(nx*(ny+1)+(ny+2)/2-1),1)=-2;
F(2*(nx*(ny+1)+(ny+2)/2),1)=-2;
F(2*(nx*(ny+1)+(ny+2)/2+1),1)=-2;
F(2*(nx*(ny+1)+(ny+2)/2+2),1)=-1;
FixDOF=1:2*(ny+1);

AllDOF=1:2*(ny+1)*(nx+1);
FreeDOF=setdiff(AllDOF,FixDOF);
% Solve linear algebra
K=sparse(K);
U(FreeDOF,:)=K(FreeDOF,FreeDOF)\F(FreeDOF,:);
U(FixDOF,:)=0;

% Mean compliance
l=full(U(:,1)'*K*U(:,1));
end

%%%%%% Sensitivity analysis %%%%%%
function [dldr]=fesens(U,top,rho)
ny=top.ny;
nx=top.nx;
pnl=top.pnl;

dldr=zeros((ny+1),(nx+1));
for iy = 1:ny
for ix = 1:nx
  n1 = (ny+1)*(ix-1)+iy;
  n2 = (ny+1)* ix +iy;
  Uel = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
  enode = [n1; n2; n2+1; n1+1] ;
  xc=[ix-1 ix   ix ix-1];
  yc=[iy-1 iy-1 iy iy];
  [DKE]=camdsen(rho(enode),pnl,xc,yc);
  for i=1:4
    dldr(enode(i)) = dldr(enode(i)) - Uel'*DKE(:,:,i)*Uel;
  end
end
end
end

%%%%%% Element stiffness matrix (Node-base method) %%%%%%
function [Ke, EE] = camd(xe,pnl,xc,yc)
[point,weight]=feglqd;
[Cm1]=fecmt1;
[Cm2]=fecmt2;
Ke=zeros(8,8);
EE=0;

for i=1:4
    [shape,dhds,dhdt]=feshape(point(i,1),point(i,2));
    [jacob]=fejacob(xc,yc,dhds,dhdt);
    [Bme]=febme(dhds,dhdt);
    
    rhoe=0;
    for j=1:4
        rhoe=rhoe+shape(j)*xe(j);
    end
    Ke = Ke + weight(i)*Bme'*((Cm1*rhoe^pnl+Cm2*(1-rhoe^pnl)))*Bme*det(jacob);
    EE = EE + rhoe*det(jacob);
end
end

%%%%%%%%%% Sensitivity: Node-base method %%%%%%%%%%
function [DK] = camdsen(xe,pnl,xc,yc)
% DK: Design sensitivity dK/dx
[point]=feglqd;
[Cm1]=fecmt1;
[Cm2]=fecmt2;
DK=zeros(8,8,4);

for i=1:4
    [shape,dhds,dhdt]=feshape(point(i,1),point(i,2));
    [jacob]=fejacob(xc,yc,dhds,dhdt);
    [Bme]=febme(dhds,dhdt);
    
    rho=0;
    for j=1:4
        rho=rho+shape(j)*xe(j);
    end
    
    for j=1:4
        C0=pnl*rho^(pnl-1)*(Cm1-Cm2);
        DK0 = Bme'*C0*Bme*det(jacob);
        DK(:,:,j) = DK(:,:,j) + shape(j)*DK0;
    end
end
end

%%%%% Gauss-Legendre quadrature %%%%%
function [point,weight]=feglqd
point=[-3^(-0.5) -3^(-0.5)
        3^(-0.5) -3^(-0.5)
        3^(-0.5)  3^(-0.5)
       -3^(-0.5)  3^(-0.5)];
weight=[1; 1; 1; 1];
end

%%%%%%%%%% Shape functions and derivatives %%%%%%%%%%
function [shape,dhds,dhdt]=feshape(s,t)
% 1st node at (-1,-1), 2nd node at ( 1,-1)
% 3rd node at ( 1, 1), 4th node at (-1, 1)
% shape functions
shape(1)=0.25*(1-s)*(1-t);
shape(2)=0.25*(1+s)*(1-t);
shape(3)=0.25*(1+s)*(1+t);
shape(4)=0.25*(1-s)*(1+t);
% derivatives
dhds(1)= -0.25*(1-t);
dhds(2)=  0.25*(1-t);
dhds(3)=  0.25*(1+t);
dhds(4)= -0.25*(1+t);
dhdt(1)= -0.25*(1-s);
dhdt(2)= -0.25*(1+s);
dhdt(3)=  0.25*(1+s);
dhdt(4)=  0.25*(1-s);
end

%%%%% Jacobian for two-dimensional mapping %%%%%
function [jacob]=fejacob(xc,yc,dhds,dhdt)
jacob=zeros(2,2);
for i=1:4
    jacob(1,1)=jacob(1,1)+dhds(i)*xc(i);
    jacob(1,2)=jacob(1,2)+dhds(i)*yc(i);
    jacob(2,1)=jacob(2,1)+dhdt(i)*xc(i);
    jacob(2,2)=jacob(2,2)+dhdt(i)*yc(i);
end
end

%%%%%%%%%% Kinematic equation between strains and displacements %%%%%%%%%%
function [Bme]=febme(dhds,dhdt)
Bme = [ dhds(1) 0       dhds(2) 0       dhds(3) 0       dhds(4) 0
        0       dhdt(1) 0       dhdt(2) 0       dhdt(3) 0       dhdt(4) 
        dhdt(1) dhds(1) dhdt(2) dhds(2) dhdt(3) dhds(3) dhdt(4) dhds(4)];
end
    
function [Cm1]=fecmt1
Emt1= 1.0;
nu1 = 0.3;
Cm1=[1  nu1  0
    nu1  1  0
    0 0 (1-nu1)/2]/(1-nu1^2)*Emt1;
end

function [Cm2]=fecmt2
Emt2= .001;
nu2 = 0.3;
Cm2=[1  nu2  0
    nu2  1  0
    0 0 (1-nu2)/2]/(1-nu2^2)*Emt2;
end
