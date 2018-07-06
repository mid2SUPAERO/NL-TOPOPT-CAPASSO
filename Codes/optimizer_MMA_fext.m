%%% Main program for the optimization of the airfoil using MMA and
%%% nonlinear mechanics
clear all
close all
clc
%% Mesh recover
load meshfina0012

% p: matrix with the coordinates of all the nodes in (x,y)
% t: connectivity matrix
p=p.';     % if the mesh comes from pdetool, we have to transpose the matrixes in order to have p:nx2 and t:ntx3
t=t.';
t=t(:,1:3);
%simpplot(p,t);              %bloc pour representer le maillage
corde = 150;
p=p*1000;           % we bring the coordinates in mm
%load resultsfrommeshfina0012P4s0_9r5
%% Extraction of geometry from the mesh
h=30;               % thickness of the plate in mm
% Coordinates of the nodes of each triangle
x1 = p(t(:,1),1);       
y1 = p(t(:,1),2);
x2 = p(t(:,2),1);
y2 = p(t(:,2),2);
x3 = p(t(:,3),1);
y3 = p(t(:,3),2);
AT = .5*(x1.*y2+y1.*x3+x2.*y3-x3.*y2-y3.*x1-x2.*y1); %area of each triangle

% Expressions used for the computation of the stiffness matrix
a2 = x1-x3;
a1 = x3-x2;
b1 = y2-y3;
a3 = x2-x1;
b3 = y1-y2;
b2 = y3-y1;
nt = length(t(:,1));    % number of triangles
% Material properties in undeformed configuration (MPa)
E0 = 15;
nu = 0.3;
lambda0 = E0*nu/(1+nu)/(1-2*nu);
mu0 = E0/2/(1+nu);
% ddl is the ntx6 matrix which in each row reports the indexes of the dofs 
%related to each triangle ordered according to the rows of the connectivity
% matrix t. The horizontal dof will be indicated as 2*node_index-1, while
% the vertical dof will be indicated as 2*node_index
ddl = [2*t(:,1)-1 2*t(:,1) 2*t(:,2)-1 2*t(:,2) 2*t(:,3)-1 2*t(:,3)];
%% Matrixes from geometry
Nux=@(b1,b2,b3,AT) [b1 0 b2 0 b3 0]/2/AT;       % shape function derivative
Nvx=@(b1,b2,b3,AT) [0 b1 0 b2 0 b3]/2/AT;
Nuy=@(a1,a2,a3,AT) [a1 0 a2 0 a3 0]/2/AT;
Nvy=@(a1,a2,a3,AT) [0 a1 0 a2 0 a3]/2/AT;
A=zeros(6,6);
AA =zeros(6,6,nt); 
B = AA;
for i=1:nt
    % A = (Nux.')*Nvy -(Nuy.')*Nvx; 
    % This is related to J=det(F) where F is the deformation tensor
    A=Nux(b1(i),b2(i),b3(i),AT(i)).'*Nvy(a1(i),a2(i),a3(i),AT(i))-Nuy(a1(i),a2(i),a3(i),AT(i)).'*Nvx(b1(i),b2(i),b3(i),AT(i));
    % B = (Nux.')*Nux + (Nuy.')*Nuy + (Nvx.')*Nvx + (Nvy.')*Nvy;
    % This is related to tr(C) where C is the right Cauchy tensor
    B(:,:,i) = (Nux(b1(i),b2(i),b3(i),AT(i)).')*Nux(b1(i),b2(i),b3(i),AT(i)) + (Nuy(a1(i),a2(i),a3(i),AT(i)).')*Nuy(a1(i),a2(i),a3(i),AT(i)) + (Nvx(b1(i),b2(i),b3(i),AT(i)).')*Nvx(b1(i),b2(i),b3(i),AT(i)) + (Nvy(a1(i),a2(i),a3(i),AT(i)).')*Nvy(a1(i),a2(i),a3(i),AT(i));
    AA(:,:,i) = A+A.';      % this is used for the successive derivatives
end
%% Boundary conditions
% We fix the nodes on the boundary of the airfoil which are in intermediate
% position
bound = boundary(p);
% figure
% plot(p(bound,1),p(bound,2));
x_bound = p(bound,1);
y_bound = p(bound,2);
%intermediate_nodes = find(p(:,1)>=0.3*corde & p(:,1)<0.75*corde);
intermediate_nodes = find(p(:,1)>=0.35*corde & p(:,1)<0.65*corde); % for meshmedia it's perfect
%intermediate_nodes = find(p(:,1)>=0.35*corde & p(:,1)<0.75*corde);
%intermediate_nodes = find(p(:,1)<=45);  % used in MBB
imposed_nodes = intersect(bound,intermediate_nodes);

figure
simpplot(p,t); hold on;
plot(p(imposed_nodes,1),p(imposed_nodes,2),'*'); hold on;

%% Nonlinear system

% ind is the index of the node where we apply the force
% It's computed as the nearest node to a point place in (0.35*corde, 0)
[dist, ind] = min(sqrt((p(:,1)-0.3*corde).^2+(p(:,2)-5).^2));
%[dist, ind] = min(sqrt((p(:,1)).^2+(p(:,2)-max(p(:,2))).^2)); %used in MBB
plot(p(ind,1),p(ind,2),'*','MarkerSize',10); hold off;
[dist, ind_borde_fuite] = min(sqrt((p(:,1)-max(p(:,1))).^2+(p(:,2)).^2));

%[dist2, ind_block] =
%min(sqrt((p(:,1)-max(p(:,1))).^2+(p(:,2)-min(p(:,2))).^2)); % used in MBB

% We block both nodal displacements of each node on the list
imposed_dofs = [2*imposed_nodes-1;2*imposed_nodes];
%imposed_dofs = [2*imposed_nodes-1; 2*ind_block];
total_dofs = [1:max(max(ddl))].';
free_dofs = setdiff(total_dofs,imposed_dofs);
%free_dofs = sort([free_dofs;2*ind-1]);
%% Mesh-independency filter preparation
rmin = 6;                       %depends on the mesh-size (it's not a contradiction with the section title)
xG = (x1+x2+x3)/3;
yG = (y1+y2+y3)/3;
i=1:nt;
[I,J]=meshgrid(i,i);
W = max(0,rmin-sqrt((xG(J)-xG(I)).^2+(yG(J)-yG(I)).^2));
W = W./repmat(sum(W,2),1,nt);
W = sparse(W);
tail_part = find(xG>0.85*corde);
%% Problem statement
% min compliance
% V<volfrac
% U(borde_fuite)<U0
U0 = -21;
volfrac = 1;
%sigmalim = 1;
sigmalim = 1;
x = volfrac*ones(nt,1);
xmin = 0.1*ones(nt,1);  % this will change at each loop iteration
xmax =ones(nt,1);
penal = 3;              % this will change at specific loop iterations
%% INITIALIZE ITERATION
%x = volfrac*ones(nt,1);
rho = x;
loop = 0;
change = 1;
m = 2;      % number of constraints
n = length(rho);    % number of variables
epsimin = 0.0000001;
eeen    = ones(n,1);
eeem    = ones(m,1);
zeron   = zeros(n,1);
zerom   = zeros(m,1);
xval    = rho;
rho = W*xval;
    figure(10)
    patch('Faces',t,'Vertices',p,'LineStyle','none','FaceVertexCData',repmat((1-rho)/(1-0.1),1,3),'FaceColor','flat');
    axis equal; drawnow;hold on;
    plot(p(imposed_nodes,1),p(imposed_nodes,2),'*'); hold on;
    plot(p(ind,1),p(ind,2),'*'); hold on;
    plot(p(ind_borde_fuite,1),p(ind_borde_fuite,2),'*');
    hold off
xold1   = xval;
xold2   = xval;
%xmin    = zeron;
xmax    = eeen;
low     = xmin;     % lower asymptote for MMA
upp     = xmax;     % upper asymptote for MMA
C       = 1000*eeem;    % proposed by Svanberg
d       = 0*eeem;
a0      = 1;
a       = zerom;
outeriter = 0;
maxoutit  = 15000;
kkttol  =0.002;
Emin = 0.001;
%%%% The iterations start:
kktnorm = kkttol+10;
% kktnorm = kkttol;
outit = 0;
cont = 0;
contmax = 200;
%load resmeshbuonaprov
%% START ITERATION
while (kktnorm > kkttol) && (outit < maxoutit) && (cont<contmax)
    outit   = outit+1;
    outeriter = outeriter+1;
    %%% Inferior limit update (expression provided by Bhattacharyya)
    
    %%% Computation of the objective function and its derivatives
    [o,do,U_15, dR_drho, K_15] = objective_fcn_fext( xval, W, penal, ddl, AT, AA, B, a1, a2, a3, b1, b2, b3, mu0, lambda0, h, imposed_dofs, ind, free_dofs);
    %%% Computation of the constraints functions and its derivatives
    %[c,~,dc,~] = constraints_fcn2( xval, W, AT, volfrac, t, p, ind_borde_fuite, U_15, dR_drho, K_15, free_dofs, U0 ); %%%with U_borde_fuite>U0
    %[ c, ceq, dc, dceq ] = constraints_fcn2( x, W, AT, volfrac,t,p,ind_borde_fuite, U_15, dR_drho, K_15, free_dofs, U0, sigmalim, nt,lambda0, mu0, a1, a2, a3, b1, b2, b3, ddl );
    [ c, ceq, dc, dceq, sigma ] = constraints_fcn2( xval, W, AT, tail_part,t,p,ind_borde_fuite, U_15, dR_drho, K_15, free_dofs, U0, sigmalim, nt,lambda0, mu0, a1, a2, a3, b1, b2, b3, ddl );
    f0val=o;
    fval=c';
    df0dx=do(:);
    dfdx=dc';
    innerit=0;
    outvector1 = [outeriter innerit f0val fval'];
    outvector2 = xval;
    %%% MMA code optimization
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,S,low,upp] = ...
        mmasub(m,n,outeriter,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,a,C,d);
    xold2 = xold1;
    xold1 = xval;
    xval  = xmma;
    %%% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Const.:%7.3f %7.3f %7.3f kktnorm.:%7.3f\n',outit,o, ...
    c(1),c(2),c(1), kktnorm);
    %%%% The residual vector of the KKT conditions is calculated:
    [residu,kktnorm,residumax] = ...
        kktcheck(m,n,xval,ymma,zmma,lam,xsi,eta,mu,zet,S, ...
        xmin,xmax,df0dx,fval,dfdx,a0,a,C,d);
    %%% PLOT DENSITIES INSIDE THE constraints_fcn
    if norm(xval-xold1)<0.0001*norm(xval)
      cont = cont+1;
    end
    rho = W*xval;
    figure(10)
    patch('Faces',t,'Vertices',p,'LineStyle','none','FaceVertexCData',repmat((1-rho)/(1-0.1),1,3),'FaceColor','flat');
    axis equal; drawnow;hold on;
    plot(p(imposed_nodes,1),p(imposed_nodes,2),'*'); hold on;
    plot(p(ind,1),p(ind,2),'*'); hold on;
    plot(p(ind_borde_fuite,1),p(ind_borde_fuite,2),'*');
    hold off
    figure(12)
    patch('Faces',t,'Vertices',p,'LineStyle','none','FaceVertexCData',sigma.','FaceColor','flat');
    axis equal; drawnow;hold on;
    plot(p(imposed_nodes,1),p(imposed_nodes,2),'*'); hold on;
    plot(p(ind,1),p(ind,2),'*'); hold on;
    plot(p(ind_borde_fuite,1),p(ind_borde_fuite,2),'*');
    colorbar
    hold off
end
%% Alternative resolution with fmincon (too slow)
% options = optimoptions('fmincon','Algorithm','interior-point','Display','iter',...
%     'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'MaxIterations',100000);
% [X] = fmincon(@(x) objective_fcn_fext( x, W, penal, ddl, AT, AA, B, a1, a2, a3, b1, b2, b3, mu0, lambda0, h, imposed_dofs, ind, free_dofs),x,[],[],[],[],xmin, xmax,@(x) constraints_fcn( x, W, AT,volfrac,t,p ),options);
%% Plot the final configuration with the relative displacement field
rho = W*xval;       % density filtering
mue = mu0 *rho.^penal;
lambdae = lambda0 *rho.^penal;
[ U_15, K_15, Fint_15] = nonlinear_tri_fext( ddl, AT, AA, B, a1, a2, a3, b1, b2, b3, mue, lambdae, h, imposed_dofs, ind );
pp=p;
for i=1:length(p(:,1))
    pp(i,1)=p(i,1)+U_15(2*i-1);
    pp(i,2)=p(i,2)+U_15(2*i);
end
figure
simpplot(pp,t); hold on;
plot(x1,y1,'*'); hold on;
plot(x2,y2,'*'); hold on;
plot(x3,y3,'*'); hold on;
plot(p(imposed_nodes,1),p(imposed_nodes,2));hold off;
[ sigmaVM, dsigmaVM_dUf, sigmax, sigmay, sigmaxy, sigmaz, dsigmax_dUf, dsigmay_dUf, dsigmaxy_dUf, dsigmaz_dUf ] = stress_fcn( U_15, ddl, free_dofs, lambda0, mu0, a1, a2, a3, b1, b2, b3, AT);
figure(13)
patch('Faces',t,'Vertices',p,'LineStyle','none','FaceVertexCData',sigmaVM,'FaceColor','flat');
axis equal; drawnow;hold on;
plot(p(imposed_nodes,1),p(imposed_nodes,2),'*'); hold on;
plot(p(ind,1),p(ind,2),'*'); hold on;
plot(p(ind_borde_fuite,1),p(ind_borde_fuite,2),'*');
colorbar
hold off