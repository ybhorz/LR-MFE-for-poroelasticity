function Err = Biot_P1(Domn,Parm,Sol,post)

%% tri-subdivision: mesh2D_h0 (h), mesh2D_h1 (h/2), mesh2D_h2 (h/4)
% % clear; close all;
% load("mesh2D_[-1,1,-1,1]\mesh2D_h0.mat");
% % p: node coordinate
% % e: boundary 
% %   -node 1; -node 2 
% %   -cumulate distance 1; -cumulate distance 2 
% %   -which boundary
% %   -?
% %   -which domain
% % t: triangle element
% %   -node 1; -node 2; node 3;
% %   -which domain
% node = p; nNode = size(node,2);
% elem = t(1:3,:); nElem = size(elem,2);
% [conn,elemCn,edge] = genConn(elem,e([1,2,5],:)); 
% nEdge = edge.nConn; nConn = conn.nConn;
% 
% EgDire = [2,3,4]; EgNeum = 1;
%% unifrom-triangulation
% domn = [0,1,0,1];
% nSub = 32;
domn = Domn.domn;
nSub = Domn.nSub;
[node,elem,edge,conn,elemCn] = genTri(domn,nSub);
nNode = size(node,2); nElem = size(elem,2); nEdge = edge.nConn; nConn = conn.nConn;

% EgDire = [2,3,4]; EgNeum = 1;
EgDire = Domn.EgDire; EgNeum = Domn.EgNeum;
%% Real solution
% lmd = 1; mu = 1; alph = 1; c0 = 0; %k = 1;
lmd = Parm.lmd; mu = Parm.mu; alph = Parm.alph; c0 = Parm.c0; k = Parm.k;

syms x y;
var = [x;y];

% ux = Fcn(var,sin(2*pi*y)*(-1+cos(2*pi*x))+1/lmd*sin(pi*x)*sin(pi*x));
% uy = Fcn(var,sin(2*pi*x)*(1-cos(2*pi*y))+1/lmd*sin(pi*x)*sin(pi*y));
eval(['ux = Fcn(var,', Sol.ux, ');']);
eval(['uy = Fcn(var,', Sol.uy, ');']);
u = [ux;uy];

% p = Fcn(var,vp*cos(pi*x)*sin(pi*y));
eval(['p = Fcn(var,', Sol.p, ');']);
z = [p.dif([1;0])*(-k);p.dif([0;1])*(-k)];

uxFun = ux.getFun;
uyFun = uy.getFun;
pFun = p.getFun;

fx = - ux.dif([2;0]) * (lmd+2*mu) - ux.dif([0;2]) * mu - uy.dif([1;1]) * (lmd+mu) + p.dif([1;0]) * alph;
fy = - uy.dif([0;2]) * (lmd+2*mu) - uy.dif([2;0]) * mu - ux.dif([1;1]) * (lmd+mu) + p.dif([0;1]) * alph;
g = p * c0 + ux.dif([1;0]) * alph + uy.dif([0;1]) * alph - p.dif([2;0]) * k - p.dif([0;2]) * k;
%% transformation
syms l m x1 x2 x3 y1 y2 y3;
refVar = [l;m];

tf = Tf;
tf.refVar = [l;m];
tf.orgVar = [x;y];
tf.parm = [x1,x2,x3;y1,y2,y3];
B = [x2-x1, x3-x1; y2-y1, y3-y1];
b = [x1;y1];
tf.toRef = B * tf.refVar + b;
tf.toOrg = B \ (tf.orgVar - b);

syms s;
tfI = Tf;
tfI.refVar = s;
tfI.orgVar = [x;y];
tfI.parm = [x1,x2;y1,y2];
tfI.toRef = [(x2-x1)/2*s+(x2+x1)/2;
            (y2-y1)/2*s+(y2+y1)/2];
%% Mesh
msh = Msh;
msh.node = node;
msh.elem = elem;
msh.edge = edge;
msh.conn = conn;
msh.tf = tf;
msh.tfI = tfI;
msh.gInt = GInt("D2P4");
msh.gIntI = GInt("D1P5");
msh.check;
%% Finite element space: P0
P0 = FE;
P0.nNode = nElem;
P0.elem = 1:nElem;
org = msh.tf.getOrg;
P0.node = zeros([2,nElem]);
for iElem = 1:nElem
    parm = node(:,elem(:,iElem));
    P0.node(:,iElem) = org([1/3;1/3],parm);
end
P0.base = Base(refVar,1,1);
P0.tf = tf; P0.baseTf = 1;
%% Finite element space: P1
P1 = FE;
P1.nNode = nNode;
P1.node = node;
P1.elem = elem;
P1.edge = edge;
P1.base = Base(refVar,[1-l-m,l,m],eye(3));
P1.tf = tf; P1.baseTf = ones(1,P1.nBase);
%% Finite element space: RT0
RT0 = FE;
RT0.nNode = conn.nConn;
RT0.node = (node(:,conn.node(1,:))+node(:,conn.node(2,:)))/2;
RT0.elem = elemCn;
RT0.edge = edge;
RT0.edge.node = 1:nEdge;
RT0.base = Base(refVar,[1,0,l;0,1,m],[0,-1,1;0,0,1;-1,0,1]');
RTtf = tf;
RTtf.refCoef = RTtf.getJdet.fun * inv(B);
RTtf.orgCoef = 1/RTtf.getJdet.fun * B;
RT0.tf = RTtf; RT0.baseTf = ones(1,RT0.nBase);
%% Finite element space
Uh = P1.^2;
Ph = P0;
Zh = RT0;
%% Varitional equation
% trls = [Uh,Zh,Ph]; tsts = trls;
iu = 1; iz = 2; ip = 3;
d0ux = [0,nan;0,nan]; d0uy = [nan,0;nan,0];
dxux = [1,nan;0,nan]; dyux = [0,nan;1,nan];
dxuy = [nan,1;nan,0]; dyuy = [nan,0;nan,1];
d0zx = [0,nan;0,nan]; d0zy = [nan,0;nan,0];
divu = [1,0;0,1]; divz = [1,0;0,1];
d0p = [0;0]; dxp = [1;0]; dyp = [0;1];

Auv(1:8) = DLF;
Auv(1) = DLF(Fcn(var,lmd+2*mu),dxux,dxux);
Auv(2) = DLF(Fcn(var,mu),dyux,dyux);
Auv(3) = DLF(Fcn(var,mu),dyux,dxuy);
Auv(4) = DLF(Fcn(var,lmd),dxux,dyuy);
Auv(5) = DLF(Fcn(var,lmd+2*mu),dyuy,dyuy);
Auv(6) = DLF(Fcn(var,mu),dxuy,dxuy);
Auv(7) = DLF(Fcn(var,mu),dxuy,dyux);
Auv(8) = DLF(Fcn(var,lmd),dyuy,dxux);

Bpv = DLF(Fcn(var,-alph),d0p,divu);

Azw(1) = DLF(Fcn(var,1/k),d0zx,d0zx);
Azw(2) = DLF(Fcn(var,1/k),d0zy,d0zy);

Bpw = DLF(Fcn(var,-1),d0p,divz);

Buq = DLF(Fcn(var,alph),divu,d0p);

Bzq = DLF(Fcn(var,1),divz,d0p);

% Apq = DLF(Fcn(var,c0),ip,ip,d0p,d0p);

Fv(1) = SLF(fx,d0ux);
Fv(2) = SLF(fy,d0uy);

Gq = SLF(g,d0p);

% Normal vector.
nm(1) = Fcn;
nm(1).var = [x;y];
nm(1).parm = [x1,x2;y1,y2];
nm(1).fun = 1/sqrt((x2-x1)^2+(y2-y1)^2) * (y2-y1);
nm(2) = Fcn;
nm(2).var = [x;y];
nm(2).parm = [x1,x2;y1,y2];
nm(2).fun = - 1/sqrt((x2-x1)^2+(y2-y1)^2) * (x2-x1);

% Neuman boundary condition.
Fv(3) = SLF;
Fv(3).load = nm(1).*ux.dif([1;0])*(lmd+2*mu) + nm(2).*ux.dif([0;1])*mu + nm(1).*uy.dif([0;1])*lmd + nm(2).*uy.dif([1;0])*mu - nm(1).*p*alph;
Fv(3).tstOrd = d0ux;
Fv(3).type = "edge";
Fv(3).domn = find(ismember(msh.edge.type,EgNeum));

Fv(4) = SLF;
Fv(4).load = nm(2).*uy.dif([0;1])*(lmd+2*mu) + nm(1).*uy.dif([1;0])*mu + nm(2).*ux.dif([1;0])*lmd + nm(1).*ux.dif([0;1])*mu - nm(2).*p*alph;
Fv(4).tstOrd = d0uy;
Fv(4).type = "edge";
Fv(4).domn = find(ismember(msh.edge.type,EgNeum));

Gz(1) = SLF;
Gz(1).load = -nm(1).*p;
Gz(1).tstOrd = d0zx;
Gz(1).type = "edge";
Gz(1).domn = find(ismember(msh.edge.type,EgNeum));

Gz(2) = SLF;
Gz(2).load = -nm(2).*p;
Gz(2).tstOrd = d0zy;
Gz(2).type = "edge";
Gz(2).domn = find(ismember(msh.edge.type,EgNeum));
%% Assmeble matrix
Stiff11 = assemble(msh,Uh,Uh,Auv,SLF.empty,"Integral",GInt('D2P1')); Stiff11 = sparse(Stiff11);
Stiff12 = sparse(Uh.nNode,Zh.nNode);
Stiff13 = assemble(msh,Ph,Uh,Bpv,SLF.empty,"Integral",GInt('D2P1')); Stiff13 = sparse(Stiff13);
[~,Load1] = assemble(msh,FE.empty,Uh,DLF.empty,Fv,"Integral",GInt('D2P4'));

Stiff21 = sparse(Zh.nNode,Uh.nNode);
Stiff22 = assemble(msh,Zh,Zh,Azw,SLF.empty,"MatrixType","sparse","Integral",GInt('D2P2')); Stiff22 = sparse(Stiff22);
Stiff23 = assemble(msh,Ph,Zh,Bpw,SLF.empty,"Integral",GInt('D2P1')); Stiff23 = sparse(Stiff23);
[~,Load2] = assemble(msh,FE.empty,Zh,DLF.empty,Gz,"Integral",GInt('D2P4'));

Stiff31 = assemble(msh,Uh,Ph,Buq,SLF.empty,"Integral",GInt('D2P1')); Stiff31 = sparse(Stiff31);
Stiff32 = assemble(msh,Zh,Ph,Bzq,SLF.empty,"Integral",GInt('D2P1')); Stiff32 = sparse(Stiff32);
Stiff33 = sparse(Ph.nNode,Ph.nNode);
[~,Load3] = assemble(msh,FE.empty,Ph,DLF.empty,Gq,"Integral",GInt('D2P4'));
%% Dirichlet boundary condition
% type 1: y = 1
% type 2: x = 1
% type 3: y = 0 
% type 4: x = 0

uxBd = uxFun; 
uyBd = uyFun; 

for iEdge = 1:Uh.nEdge 
    for iNode = 1:2
        I = Uh.edge.node(iNode,iEdge);
        coord = Uh.node(:,I);
        if ismember(Uh.edge.type(iEdge),EgDire)
            Stiff11(I,:) = 0;
            Stiff11(I,I) = 1;
            Stiff13(I,:) = 0;
            Load1(I) = uxBd(coord);
        end
    end
    for iNode = 3:4
        I = Uh.edge.node(iNode,iEdge);
        coord = Uh.node(:,I);
        if ismember(Uh.edge.type(iEdge),EgDire)
            Stiff11(I,:) = 0;
            Stiff11(I,I) = 1;
            Stiff13(I,:) = 0;
            Load1(I) = uyBd(coord);
        end
    end
end

znBd = -nm(1).*p.dif(dxp) - nm(2).*p.dif(dyp);
znBd = znBd.getFun;
orgI = msh.tfI.getOrg;
JnormI = msh.tfI.getJnorm.getFun;

for iEdge = 1:Zh.nEdge
    EgParm = msh.node(:,msh.edge.node(:,iEdge));
    I = Zh.edge.node(iEdge);
    if ismember(Zh.edge.type(iEdge),EgDire)
        Stiff22(I,:) = 0;
        Stiff22(I,I) = 1;
        Stiff23(I,:) = 0;
        Load2(I) = msh.gIntI.getVal(@(t) znBd(orgI(t,EgParm),EgParm).* JnormI(t,EgParm));
    end
end

% Fix p in one point
if isempty(EgNeum)
    I = 1;
    coord = Ph.node(:,I);
    Stiff31(I,:) = 0;
    Stiff32(I,:) = 0;
    Stiff33(I,I) = 1;
    Load3(I) = pFun(coord);
end

Stiff = [Stiff11,Stiff12,Stiff13;
        Stiff21,Stiff22,Stiff23;
        Stiff31,Stiff32,Stiff33];
Load = [Load1;Load2;Load3];
%% Solve equation
xh = full(Stiff\Load);

nTrlNodes = cumsum([0,Uh.nNode,Zh.nNode,Ph.nNode]);
uh = xh(nTrlNodes(iu)+1:nTrlNodes(iu+1));
zh = xh(nTrlNodes(iz)+1:nTrlNodes(iz+1));
ph = xh(nTrlNodes(ip)+1:nTrlNodes(ip+1));
%% Error norms
eu_h0 = eNorm(msh,Uh,uh,u,{d0ux,d0uy}); 
eu_h1 = eNorm(msh,Uh,uh,u,{dxux,dyux,dxuy,dyuy});
ez_h0 = eNorm(msh,Zh,zh,z,{d0zx,d0zy});
ez_div = eNorm(msh,Zh,zh,z,divz);
ep_h0 = eNorm(msh,Ph,ph,p,d0p);
[~,Php] = assemble(msh,FE.empty,P0,DLF.empty,SLF(p*2,1,d0p,'refElem'),"Integral",GInt("D2P4"));
ePhp_h0 = eNorm(msh,Ph,Php-ph,Fcn(var,0),d0p);

fprintf('h: %f e(u)_h0: %f e(u)_h1: %f e(z)_h0: %f e(z)_div: %f e(p)_h0: %f e(Php)_h0: %f\n',msh.h,eu_h0,eu_h1,ez_h0,ez_div,ep_h0,ePhp_h0);
%%
Err = table(eu_h0,eu_h1,ez_h0,ez_div,ep_h0,ePhp_h0);
end