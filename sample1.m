% Topology Optimization
function [l,v]=sample1

% Problem Settings
top.nx = 120;   % The number of elements in the horizontal axis
top.ny = 60;   % The number of elements in the vertical axis
top.vol  = .5;  % Volume constraint
top.pnl  = 3;   % Penalization parameter
top.rfil = 1.5; % Filtering distance
% Parameters for optimization
opt.sw1=0;  % OC
%opt.sw1=1; % SLP
opt.convprm1=0.01;
opt.convprm2=20;

% Initialization
rho=top.vol*ones(top.ny,top.nx);
itr = 0;
minl=1.e30;
convflag=0;

% Display initial configuration
colormap(gray);imagesc(1-rho);caxis([0 1]);
axis equal;axis off;drawnow;

% Optimization loop begins here
while convflag <= opt.convprm2
    itr = itr + 1;
    [U, l]=fea(rho,top);            % Finite element analysis
    [dldr]=fesens(U,rho,top);       % Sensitivity calculation
    [dldr]=Filter(top,rho,dldr);    % Sensitivity filtering

    if opt.sw1==0 % Design variable update
        [rhonew] = OC(top,rho,dldr);
    else
        [rhonew] = SLP(top,rho,dldr);
    end
    
    % Convergence check
    change = max(max(abs(rhonew-rho)));
    if minl>l
        minl=l;
        convflag=0;
    else
        convflag=convflag+1;
    end
    if change < opt.convprm1
        convflag=opt.convprm2+1;
    end
    rho=rhonew;
    
    v=mean(mean(rhonew));
    % Show optimization process data
    fprintf('It.:%4i Cmp.:%8.3f Chng.:%7.3f Vol.:%6.3f\n', itr, l, change, v);
    imagesc(1-rho);
    axis equal;axis off;drawnow;    % Show density distribution

end
end

%%%%%%%%%% Optimality Criteria %%%%%%%%%
function [rhonew]=OC(top,rho,dldr)
vol=top.vol;
rho=reshape(rho,top.nx*top.ny,1);
dldr=reshape(dldr,top.nx*top.ny,1);
lambda1 = 0; lambda2 = 1e4;
mvlmt = 0.15; %Move-limit
eta=0.5;  %Damping factor
while (lambda2-lambda1)/(lambda1+lambda2) > 1e-3
    lmid = (lambda2+lambda1)*0.5;
    B=rho.*(-dldr./lmid).^eta;
    rhonew = max(max(min(min(B,rho+mvlmt),1),rho-mvlmt),0);
    if (mean(rhonew(:))-vol>0)
        lambda1=lmid; 
    else
        lambda2=lmid;
    end
end
rhonew = reshape(rhonew,top.ny,top.nx);
end

%%%%%%%%%% Sequential Linear Programming %%%%%%%%%
function [rhonew]=SLP(top,rho,dldr)
nx=top.nx;
ny=top.ny;
vol=top.vol;
nvar=nx*ny; %The number of variables
mvlmt = 0.1; %Move-limit
LB = 0.;    %Lower bound of variables
UB = 1.;       %Upper bound of variables
rho=reshape(rho,nvar,1);

%Side constraints
LBm=max(LB,rho*(1-mvlmt));
UBm=min(UB,rho*(1+mvlmt));

%Volume constraint
b(1)=1;
A=ones(1,nvar)/(vol*nvar);

%Solve linear programming problem
rhonew = reshape(linprog(dldr,A,b,[],[],LBm,UBm),ny,nx);
end

%%%%%%%%%% Sensitivity filtering %%%%%%%%%%%
function [dlf]=Filter(top,rho,dldr)
nx=top.nx;
ny=top.ny;
rfil=top.rfil;
dlf=zeros(ny,nx);
for ix = 1:nx
for iy = 1:ny
    sum=0.0;
    for i = max(ix-ceil(rfil),1):min(ix+ceil(rfil), nx)
    for j = max(iy-ceil(rfil),1):min(iy+ceil(rfil), ny)
        fac = rfil-sqrt((ix-i)^2+(iy-j)^2);
        sum = sum+max(0,fac);
        dlf(iy,ix)=dlf(iy,ix)+max(0,fac)*rho(j,i)*dldr(j,i);
    end
    end
    dlf(iy,ix)=dlf(iy,ix)/(max(rho(iy,ix)*sum,1e-10));
end
end
end

%%%%%%%%%% Finite element anlysis %%%%%%%%%%%%
function [U,l]=fea(rho,top)
pnl=top.pnl;
nx=top.nx;
ny=top.ny;
K=zeros(2*(nx+1)*(ny+1),2*(nx+1)*(ny+1));
U=sparse(2*(ny+1)*(nx+1),1);
for iy=1:ny
for ix=1:nx
    n1=(ny+1)*(ix-1)+iy;
    n2=(ny+1)*ix+iy;
    elem=[2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
    xc=[ix-1 ix   ix ix-1];
    yc=[iy-1 iy-1 iy iy];
    [Ke]=feelk(rho(iy,ix),pnl,xc,yc);
    K(elem,elem)=K(elem,elem) + Ke;
end
end

% Boundary conditions
%Half Messerschmitt-Bolkow-Blohm beams
F=sparse(2,1,-1,2*(ny+1)*(nx+1),1);
FixDOF = union(1:2:2*(ny+1),2*(nx+1)*(ny+1));

%Short cantilever
% F=sparse(2*(nx*(ny+1)+(ny+2)/2),1,-1,2*(ny+1)*(nx+1),1);
% FixDOF=1:2*(ny+1);

%Bridge
% F=sparse(2*(nx/2+1)*(ny+1),1,1,2*(ny+1)*(nx+1),1);
% FixDOF=[2*ny+1 2*(ny+1) 2*(nx+1)*(ny+1)];

AllDOF=1:2*(ny+1)*(nx+1);
FreeDOF=setdiff(AllDOF,FixDOF);

% Solve linear algebra
K=sparse(K);
U(FreeDOF,:)=K(FreeDOF,FreeDOF)\F(FreeDOF,:);
U(FixDOF,:)= 0;

% Mean compliance
l = full(U'*K*U);
end

%%%%%%%%%% Element stiffness matrix %%%%%%%
function [Ke]=feelk(rho,pnl,xc,yc)
Emin=1e-3;
E0=1.;
Emt0= (E0-Emin)*rho^pnl+Emin;
Cm=fecmt*Emt0;
[point,weight]=feglqd;
Ke=zeros(8,8);
for i=1:4
    [dhds,dhdt]=feshape(point(i,1),point(i,2));
    [jcb]=fejcb(xc,yc,dhds,dhdt);
    [Bme]=febme(dhds,dhdt);
    Ke = Ke + weight(i)*Bme'*Cm*Bme*det(jcb);
end
end

%%%%% Gauss-Legendre quadrature points and weights %%%%%
function [point,weight]=feglqd
point=[-3^(-0.5) -3^(-0.5)
        3^(-0.5) -3^(-0.5)
        3^(-0.5)  3^(-0.5)
       -3^(-0.5)  3^(-0.5)];
weight=[1; 1; 1; 1];
end

%%%%% Shape functions, and their derivatives for isoparametric four-node %%%%%
function [dhds,dhdt,shape]=feshape(s,t)
% 1st node at (-1,-1), 2nd node at ( 1,-1)
% 3rd node at ( 1, 1), 4th node at (-1, 1)
% Shape functions
shape(1)=0.25*(1-s)*(1-t);
shape(2)=0.25*(1+s)*(1-t);
shape(3)=0.25*(1+s)*(1+t);
shape(4)=0.25*(1-s)*(1+t);
% Derivatives
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
function [jcb]=fejcb(xc,yc,dhds,dhdt)
jcb=zeros(2,2);
for i=1:4
    jcb(1,1)=jcb(1,1)+dhds(i)*xc(i);
    jcb(1,2)=jcb(1,2)+dhds(i)*yc(i);
    jcb(2,1)=jcb(2,1)+dhdt(i)*xc(i);
    jcb(2,2)=jcb(2,2)+dhdt(i)*yc(i);
end
end

%%%%% Kinematic equation between strains and displacements %%%%%
function [Bme]=febme(dhds,dhdt)
Bme = [ dhds(1) 0       dhds(2) 0       dhds(3) 0       dhds(4) 0
        0       dhdt(1) 0       dhdt(2) 0       dhdt(3) 0       dhdt(4) 
        dhdt(1) dhds(1) dhdt(2) dhds(2) dhdt(3) dhds(3) dhdt(4) dhds(4)];
end
    
%%%%% Sensitivity calculation %%%%%
function [dldr]=fesens(U,rho,top)
ny=top.ny;
nx=top.nx;
pnl=top.pnl;
dldr=zeros(ny,nx);
for iy = 1:ny
for ix = 1:nx
    n1 = (ny+1)*(ix-1)+iy;
    n2 = (ny+1)* ix +iy;
    Uel = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
    xc=[ix-1 ix   ix ix-1];
    yc=[iy-1 iy-1 iy iy];
    [DKE] = dkedr(rho(iy,ix),pnl,xc,yc);
    dldr(iy,ix) = -Uel'*DKE*Uel;
end
end
end

%%%%% Sensitivity of element stiffness matrix %%%%%
function [DKE]=dkedr(rho,pnl,xc,yc)
Emin=1e-3;
E0=1.;
dEdx= (E0-Emin)*pnl*rho^(pnl-1);
Cm=fecmt*dEdx;
[point,weight]=feglqd;
DKE=zeros(8,8);
for i=1:4
    [dhds,dhdt]=feshape(point(i,1),point(i,2));
    [jcb]=fejcb(xc,yc,dhds,dhdt);
    [Bme]=febme(dhds,dhdt);
    DKE = DKE + weight(i)*Bme'*Cm*Bme*det(jcb);
end
end

%%%%% Elastic coefficient matrix %%%%%
function [Cm]=fecmt
nu = 0.3;
Cm=[1   nu  0
    nu  1   0
    0   0   (1-nu)/2] / (1-nu^2);
end
