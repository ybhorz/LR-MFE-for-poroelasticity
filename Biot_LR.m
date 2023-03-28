function Err = Biot_LR(Domn,Parm,Sol,post)

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
%% Finite element space: P1D
P1D = FE;
P1D.nNode = 3*nElem;
P1D.elem = reshape(1:P1D.nNode,[3,nElem]);
P1D.node = zeros(2,P1D.nNode);
for iElem = 1:P1D.nElem
    for iNode = 1:3
        P1D.node(:,P1D.elem(iNode,iElem)) = node(:,elem(iNode,iElem));
    end
end

P1D.edge = edge;
for iEdge = 1:nEdge
    iElem = edge.elem(iEdge);
    switch find(elemCn(:,iElem)==iEdge)
        case 1
            P1D.edge.node(:,iEdge) = P1D.elem([1,2],iElem);
        case 2
            P1D.edge.node(:,iEdge) = P1D.elem([2,3],iElem);
        case 3
            P1D.edge.node(:,iEdge) = P1D.elem([3,1],iElem);
    end
end

% !Note: Elements with one edge vertex are also included.
for iEdge = 1:nEdge
    for iElem = setdiff(1:nElem,P1D.edge.elem)
        for iNode = 1:3
            if ismember(elem(iNode,iElem),edge.node(:,iEdge))
                P1D.edge.node = [P1D.edge.node, [P1D.elem(iNode,iElem);0]];
                P1D.edge.elem = [P1D.edge.elem, iElem];
                P1D.edge.type = [P1D.edge.type, edge.type(iEdge)];
                break;
            end
        end
    end
end

P1D.base = Base(refVar,[1-l-m,l,m],eye(3));

P1D.tf = tf; P1D.baseTf = ones(1,3);
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
%% Finite element space RT1
RT1 = FE;
RT1.nNode = 2*nConn + 2*nElem;

RT1.node = zeros(2,RT1.nNode);
RT1.node(:,1:nConn) = (1/2+1/(2*sqrt(3)))*node(:,conn.node(1,:))...
    + (1/2-1/(2*sqrt(3)))*node(:,conn.node(2,:));
RT1.node(:,nConn+(1:nConn)) = (1/2-1/(2*sqrt(3)))*node(:,conn.node(1,:))...
    + (1/2+1/(2*sqrt(3)))*node(:,conn.node(2,:));

org = msh.tf.getOrg;
for iElem = 1:nElem
    parm = node(:,elem(:,iElem));
    RT1.node(:,2*nConn+iElem) = org([1/3;1/3],parm);
    RT1.node(:,2*nConn+nElem+iElem) = org([1/3;1/3],parm);
end

RT1.elem = zeros(8,nElem);
for iElem = 1:nElem
    for iConn = 1:3
        if elemCn(iConn,iElem) > 0
            RT1.elem(iConn,iElem) = elemCn(iConn,iElem);
            RT1.elem(3+iConn,iElem) = elemCn(iConn,iElem) + nConn;
        else
            %! Note
            RT1.elem(iConn,iElem) = elemCn(iConn,iElem) - nConn;
            RT1.elem(3+iConn,iElem) = elemCn(iConn,iElem);
        end
    end
end
RT1.elem(7,:) = 2*nConn + (1:nElem);
RT1.elem(8,:) = 2*nConn + nElem + (1:nElem);

RT1.edge = edge;
RT1.edge.node = [1:nEdge;nConn+(1:nEdge)];

baseFs = [1,l,m,l*(l+m),0,0,0,0;0,0,0,0,1,m,l,m*(l+m)]; % Finite fucntion space
FsFun = matlabFunction(baseFs,'Vars',refVar);
DoFs = zeros(8); % Degree of freedom

DoFs(1,:) = [0,-1] * FsFun(1/2-1/(2*sqrt(3)),0);
DoFs(2,:) = [1/sqrt(2),1/sqrt(2)] * FsFun(1/2+1/(2*sqrt(3)),1/2-1/(2*sqrt(3)));
DoFs(3,:) = [-1,0] * FsFun(0,1/2+1/(2*sqrt(3)));
DoFs(4,:) = [0,-1] * FsFun(1/2+1/(2*sqrt(3)),0);
DoFs(5,:) = [1/sqrt(2),1/sqrt(2)] * FsFun(1/2-1/(2*sqrt(3)),1/2+1/(2*sqrt(3)));
DoFs(6,:) = [-1,0] * FsFun(0,1/2-1/(2*sqrt(3)));
DoFs(7,:) = [1,0] * FsFun(1/3,1/3);
DoFs(8,:) = [0,1] * FsFun(1/3,1/3);

baseCoef = DoFs\diag([1,1/sqrt(2),1,1,1/sqrt(2),1,1,1]); %! coefficents of base.

RT1.base = Base(refVar,baseFs,baseCoef);

RT1Tf = tf;
RT1Tf.refCoef = RT1Tf.getJdet.fun * inv(B);
RT1Tf.orgCoef = 1/RT1Tf.getJdet.fun * B;
RT1.tf = RT1Tf; RT1.baseTf = ones(1,8);
%% Finite element space
Uh1 = P1.^2;
UhR = RT0;
Ph = P0;
Zh = RT0;
%% Varitional equation
% trls = [Uh1,UhR,Zh,Ph]; tsts = trls;
iu1 = 1; iuR = 2; iz = 3; ip = 4;
d0ux = [0,nan;0,nan]; d0uy = [nan,0;nan,0];
dxux = [1,nan;0,nan]; dyux = [0,nan;1,nan];
dxuy = [nan,1;nan,0]; dyuy = [nan,0;nan,1];
d0zx = [0,nan;0,nan]; d0zy = [nan,0;nan,0];
divu = [1,0;0,1]; divz = [1,0;0,1];
d0p = [0;0]; dxp = [1;0]; dyp = [0;1];

Auv1(1:8) = DLF;
Auv1(1) = DLF(Fcn(var,lmd+2*mu),dxux,dxux);
Auv1(2) = DLF(Fcn(var,mu),dyux,dyux);
Auv1(3) = DLF(Fcn(var,mu),dyux,dxuy);
Auv1(4) = DLF(Fcn(var,lmd),dxux,dyuy);
Auv1(5) = DLF(Fcn(var,lmd+2*mu),dyuy,dyuy);
Auv1(6) = DLF(Fcn(var,mu),dxuy,dxuy);
Auv1(7) = DLF(Fcn(var,mu),dxuy,dyux);
Auv1(8) = DLF(Fcn(var,lmd),dyuy,dxux);

Bpv1 = DLF(Fcn(var,-alph),d0p,divu);

sgm = 1; % Parameter
AuvR(1) = DLF(Fcn(var,sgm/(msh.h^2)),d0ux,d0ux);
AuvR(2) = DLF(Fcn(var,sgm/(msh.h^2)),d0uy,d0uy);
AuvR(3) = DLF(Fcn(var,lmd),divu,divu);

Au1vR = DLF(Fcn(var,lmd),divu,divu);
AuRv1 = DLF(Fcn(var,lmd),divu,divu);

BpvR = DLF(Fcn(var,-alph),d0p,divu);

Azw(1) = DLF(Fcn(var,1/k),d0zx,d0zx);
Azw(2) = DLF(Fcn(var,1/k),d0zy,d0zy);

Bpw = DLF(Fcn(var,-1),d0p,divz);

Buq = DLF(Fcn(var,alph),divu,d0p);
Buq1 = DLF(Fcn(var,alph),divu,d0p);
BuqR = DLF(Fcn(var,alph),divu,d0p);

Bzq = DLF(Fcn(var,1),divz,d0p);

% Apq = DLF(Fcn(var,c0),ip,ip,d0p,d0p);

Fv1(1) = SLF(fx,d0ux);
Fv1(2) = SLF(fy,d0uy);
FvR(1) = SLF(fx,d0ux);
FvR(2) = SLF(fy,d0uy);

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
Fv1(3) = SLF;
Fv1(3).load = nm(1).*ux.dif([1;0])*(lmd+2*mu) + nm(2).*ux.dif([0;1])*mu + nm(1).*uy.dif([0;1])*lmd + nm(2).*uy.dif([1;0])*mu - nm(1).*p*alph;
Fv1(3).tstOrd = d0ux;
Fv1(3).type = "edge";
Fv1(3).domn = find(ismember(msh.edge.type,EgNeum));

Fv1(4) = SLF;
Fv1(4).load = nm(2).*uy.dif([0;1])*(lmd+2*mu) + nm(1).*uy.dif([1;0])*mu + nm(2).*ux.dif([1;0])*lmd + nm(1).*ux.dif([0;1])*mu - nm(2).*p*alph;
Fv1(4).tstOrd = d0uy;
Fv1(4).type = "edge";
Fv1(4).domn = find(ismember(msh.edge.type,EgNeum));

FvR(3) = SLF;
FvR(3).load = nm(1).*ux.dif([1;0])*(lmd+2*mu) + nm(2).*ux.dif([0;1])*mu + nm(1).*uy.dif([0;1])*lmd + nm(2).*uy.dif([1;0])*mu - nm(1).*p*alph;
FvR(3).tstOrd = d0ux;
FvR(3).type = "edge";
FvR(3).domn = find(ismember(msh.edge.type,EgNeum));

FvR(4) = SLF;
FvR(4).load = nm(2).*uy.dif([0;1])*(lmd+2*mu) + nm(1).*uy.dif([1;0])*mu + nm(2).*ux.dif([1;0])*lmd + nm(1).*ux.dif([0;1])*mu - nm(2).*p*alph;
FvR(4).tstOrd = d0uy;
FvR(4).type = "edge";
FvR(4).domn = find(ismember(msh.edge.type,EgNeum));

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
Stiff11 = assemble(msh,Uh1,Uh1,Auv1,SLF.empty,"Integral",GInt('D2P1')); Stiff11 = sparse(Stiff11);
Stiff12 = assemble(msh,UhR,Uh1,AuRv1,SLF.empty,"Integral",GInt('D2P1'),"MatrixType","sparse");
Stiff13 = sparse(Uh1.nNode,Zh.nNode);
Stiff14 = assemble(msh,Ph,Uh1,Bpv1,SLF.empty,"Integral",GInt('D2P1')); Stiff14 = sparse(Stiff14);
[~,Load1] = assemble(msh,FE.empty,Uh1,DLF.empty,Fv1,"Integral",GInt('D2P4'));

Stiff21 = assemble(msh,Uh1,UhR,Au1vR,SLF.empty,"Integral",GInt('D2P1'),"MatrixType","sparse");
Stiff22 = assemble(msh,UhR,UhR,AuvR,SLF.empty,"Integral",GInt('D2P2'),"MatrixType","sparse");
Stiff23 = sparse(UhR.nNode,Zh.nNode);
Stiff24 = assemble(msh,Ph,UhR,BpvR,SLF.empty,"Integral",GInt('D2P1')); Stiff24 = sparse(Stiff24);
[~,Load2] = assemble(msh,FE.empty,UhR,DLF.empty,FvR,"Integral",GInt('D2P4'));

Stiff31 = sparse(Zh.nNode,Uh1.nNode);
Stiff32 = sparse(Zh.nNode,UhR.nNode);
Stiff33 = assemble(msh,Zh,Zh,Azw,SLF.empty,"MatrixType","sparse","Integral",GInt('D2P2')); Stiff33 = sparse(Stiff33);
Stiff34 = assemble(msh,Ph,Zh,Bpw,SLF.empty,"Integral",GInt('D2P1')); Stiff34 = sparse(Stiff34);
[~,Load3] = assemble(msh,FE.empty,Zh,DLF.empty,Gz,"Integral",GInt('D2P4'));

Stiff41 = assemble(msh,Uh1,Ph,Buq1,SLF.empty,"Integral",GInt('D2P1')); Stiff41 = sparse(Stiff41);
Stiff42 = assemble(msh,UhR,Ph,BuqR,SLF.empty,"Integral",GInt('D2P1')); Stiff42 = sparse(Stiff42);
Stiff43 = assemble(msh,Zh,Ph,Bzq,SLF.empty,"Integral",GInt('D2P1')); Stiff43 = sparse(Stiff43);
Stiff44 = sparse(Ph.nNode,Ph.nNode);
[~,Load4] = assemble(msh,FE.empty,Ph,DLF.empty,Gq,"Integral",GInt('D2P4'));
%% Dirichlet boundary condition
% type 1: y = 1
% type 2: x = 1
% type 3: y = 0 
% type 4: x = 0

uxBd = uxFun; 
uyBd = uyFun; 

for iEdge = 1:Uh1.nEdge 
    for iNode = 1:2
        I = Uh1.edge.node(iNode,iEdge);
        coord = Uh1.node(:,I);
        if ismember(Uh1.edge.type(iEdge),EgDire)
            Stiff11(I,:) = 0;
            Stiff11(I,I) = 1;
            Stiff12(I,:) = 0;
            Stiff14(I,:) = 0;
            Load1(I) = uxBd(coord);
        end
    end
    for iNode = 3:4
        I = Uh1.edge.node(iNode,iEdge);
        coord = Uh1.node(:,I);
        if ismember(Uh1.edge.type(iEdge),EgDire)
            Stiff11(I,:) = 0;
            Stiff11(I,I) = 1;
            Stiff12(I,:) = 0;
            Stiff14(I,:) = 0;
            Load1(I) = uyBd(coord);
        end
    end
end

for iEdge = 1:UhR.nEdge
    I = UhR.edge.node(iEdge);
    Stiff21(I,:) = 0;
    Stiff22(I,:) = 0;
    Stiff22(I,I) = 1;
    Stiff24(I,:) = 0;
    Load2(I) = 0;
end

znBd = -nm(1).*p.dif(dxp) - nm(2).*p.dif(dyp);
znBd = znBd.getFun;
orgI = msh.tfI.getOrg;
JnormI = msh.tfI.getJnorm.getFun;

for iEdge = 1:Zh.nEdge
    EgParm = msh.node(:,msh.edge.node(:,iEdge));
    I = Zh.edge.node(iEdge);
    if ismember(Zh.edge.type(iEdge),EgDire)
        Stiff33(I,:) = 0;
        Stiff33(I,I) = 1;
        Stiff34(I,:) = 0;
        Load3(I) = msh.gIntI.getVal(@(t) znBd(orgI(t,EgParm),EgParm).* JnormI(t,EgParm));
    end
end

% Fix p in one point
if isempty(EgNeum)
    I = 1;
    coord = Ph.node(:,I);
    Stiff41(I,:) = 0;
    Stiff42(I,:) = 0;
    Stiff43(I,:) = 0;
    Stiff44(I,I) = 1;
    Load4(I) = pFun(coord);
end

Stiff = [Stiff11,Stiff12,Stiff13,Stiff14;
        Stiff21,Stiff22,Stiff23,Stiff24;
        Stiff31,Stiff32,Stiff33,Stiff34;
        Stiff41,Stiff42,Stiff43,Stiff44];
Load = [Load1;Load2;Load3;Load4];
%% Solve equation
xh = full(Stiff\Load);

nTrlNodes = cumsum([0,Uh1.nNode,UhR.nNode,Zh.nNode,Ph.nNode]);
uh1 = xh(nTrlNodes(iu1)+1:nTrlNodes(iu1+1));
uhR = xh(nTrlNodes(iuR)+1:nTrlNodes(iuR+1));
zh = xh(nTrlNodes(iz)+1:nTrlNodes(iz+1));
ph = xh(nTrlNodes(ip)+1:nTrlNodes(ip+1));
%% Error norms
% Uh = Uh1 + UhR; uh = [uh1;uhR];
Uh = Uh1; uh = uh1;
eu_h0 = eNorm(msh,Uh,uh,u,{d0ux,d0uy}); 
eu_h1 = eNorm(msh,Uh,uh,u,{dxux,dyux,dxuy,dyuy});
ez_h0 = eNorm(msh,Zh,zh,z,{d0zx,d0zy});
ez_div = eNorm(msh,Zh,zh,z,divz);
ep_h0 = eNorm(msh,Ph,ph,p,d0p);
[~,Php] = assemble(msh,FE.empty,P0,DLF.empty,SLF(p*2,1,d0p,'refElem'),"Integral",GInt("D2P4"));
ePhp_h0 = eNorm(msh,Ph,Php-ph,Fcn(var,0),d0p);
%% Post-processing: 1
if post == 1
    Zh = RT1; Ph = P1D; Uh = Uh1 + UhR; uh = [uh1;uhR];
    iz = 1; ip = 2;
    Stiff11 = assemble(msh,Zh,Zh,Azw,SLF.empty,"MatrixType","sparse"); Stiff11 = sparse(Stiff11);
    Stiff12 = assemble(msh,Ph,Zh,Bpw,SLF.empty,"MatrixType","sparse"); Stiff12 = sparse(Stiff12);
    Stiff20 = assemble(msh,Uh,Ph,Buq,SLF.empty,"MatrixType","sparse"); Stiff20 = sparse(Stiff20);
    Stiff21 = assemble(msh,Zh,Ph,Bzq,SLF.empty,"MatrixType","sparse"); Stiff21 = sparse(Stiff21);
    Stiff22 = sparse(Ph.nNode,Ph.nNode);
    [~,Load1] = assemble(msh,FE.empty,Zh,DLF.empty,Gz,"Integral",GInt('D2P4'));
    [~,Load2] = assemble(msh,FE.empty,Ph,DLF.empty,Gq);
    Load2 = Load2 - Stiff20*uh;
    
    for iEdge = 1:Zh.nEdge
        for iNode = 1:2
        EgParm = msh.node(:,msh.edge.node(:,iEdge));
            I = Zh.edge.node(iNode,iEdge);
            if ismember(Zh.edge.type(iEdge),EgDire)
                Stiff11(I,:) = 0;
                Stiff11(I,I) = 1;
                Stiff12(I,:) = 0;
                Load1(I) = msh.gIntI.getVal(@(t) znBd(orgI(t,EgParm),EgParm).* JnormI(t,EgParm));
            end
        end
    end
    
    % Fix p in one point
    if isempty(EgNeum)
        I = 1;
        coord = Ph.node(:,I);
        Stiff21(I,:) = 0;
        Stiff22(I,I) = 1;
        Load2(I) = pFun(coord);
    end
    
    Stiff = [Stiff11,Stiff12;Stiff21,Stiff22];
    Load = [Load1;Load2];
    xh = Stiff\Load;
    
    nTrlNodes = cumsum([0,Zh.nNode,Ph.nNode]);
    zh = xh(nTrlNodes(iz)+1:nTrlNodes(iz+1));
    ph = xh(nTrlNodes(ip)+1:nTrlNodes(ip+1));
    
    ez_h0 = eNorm(msh,Zh,zh,z,{d0zx,d0zy});
    ez_div = eNorm(msh,Zh,zh,z,divz);
    ep_h0 = eNorm(msh,Ph,ph,p,d0p);
    
    ePhp_h0 = 0;
end
%% Post-processsing: 2
if post == 2
    Ph_ = P1D; Ph = P0; Zh = RT0;
    dxp = [1;0]; dyp = [0;1];
    Apq(1) = DLF(Fcn(var,1),dxp,dxp);
    Apq(2) = DLF(Fcn(var,1),dyp,dyp);
    Bzq(1) = DLF(Fcn(var,1/k),d0zx,dxp);
    Bzq(2) = DLF(Fcn(var,1/k),d0zy,dyp);
    Cpq = DLF(Fcn(var,1),d0p,d0p);
    
    Stiff11 = assemble(msh,Ph_,Ph_,Apq,SLF.empty,"Integral",GInt("D2P1"),"MatrixType","sparse"); Stiff11 = sparse(Stiff11);
    Stiff01 = assemble(msh,Zh,Ph_,Bzq,SLF.empty,"Integral",GInt("D2P1"),"MatrixType","sparse"); Stiff01 = sparse(Stiff01);
    Load1 = - Stiff01 * zh;
    Stiff21 = assemble(msh,Ph_,Ph,Cpq,SLF.empty,"Integral",GInt("D2P1"),"MatrixType","sparse"); Stiff21 = sparse(Stiff21);
    Stiff02 = assemble(msh,Ph,Ph,Cpq,SLF.empty,"Integral",GInt("D2P1")); Stiff02 = sparse(Stiff02);
    Load2 = Stiff02 * ph;
    
    Stiff = [Stiff11;Stiff21];
    Load = [Load1;Load2];
    ph = Stiff \ Load;
    
    ep_h0 = eNorm(msh,Ph_,ph,p,d0p);
    ePhp_h0 = 0;
end
%% Plot: p
% figure; view(3);
% for iElem = 1:msh.nElem
%     phVal = repmat(ph(Ph.elem(:,iElem)),[3,1]);
%     patch('Faces',1:3,'Vertices',[msh.node(:,msh.elem(:,iElem))',phVal],'FaceVertexCData',phVal,'FaceColor','interp');
% end
% title('Approximation: p');
% 
% figure; view(3);
% for iElem = 1:nElem
%     pVal = repmat(pFun(Ph.node(:,Ph.elem(:,iElem))),[3,1]);
%     patch('Faces',1:3,'Vertices',[msh.node(:,msh.elem(:,iElem))',pVal],'FaceVertexCData',pVal,'FaceColor','interp');
% end
% title('Real solution: p');
%% Plot: u
% nUhxNodes = 1:msh.nNode;
% nUhyNodes = msh.nNode + (1:msh.nNode);
% 
% [X,Y] = meshgrid(domn(1):msh.h:domn(2),domn(3):msh.h:domn(4));
% nX = size(X,2);
% nY = size(Y,1);
% 
% uhx = uh1(nUhxNodes);
% uhy = uh1(nUhyNodes);
% uhxNode = Uh1.node(:,nUhxNodes);
% uhyNode = Uh1.node(:,nUhyNodes);
% uhxVal = griddata(uhxNode(1,:)',uhxNode(2,:)',uhx,X,Y,'cubic');
% uhyVal = griddata(uhyNode(1,:)',uhyNode(2,:)',uhy,X,Y,'cubic');
% 
% figure;
% for i = 1:nX-1
%     for j = 1:nY-1
%         vertX = [X(i,j)+uhxVal(i,j); X(i+1,j)+uhxVal(i+1,j); X(i+1,j+1)+uhxVal(i+1,j+1); X(i,j+1)+uhxVal(i,j+1)];
%         vertY = [Y(i,j)+uhyVal(i,j); Y(i+1,j)+uhyVal(i+1,j); Y(i+1,j+1)+uhyVal(i+1,j+1); Y(i,j+1)+uhyVal(i,j+1)];
%         patch('Faces',1:4,'Vertices',[vertX,vertY],'EdgeColor','black','FaceColor','none');
%     end
% end
% title('Approximation: u');
% 
% uxVal = zeros(size(X));
% uyVal = zeros(size(X));
% for i = 1:size(X,1)*size(X,2)
%     uxVal(i) = uxFun([X(i);Y(i)]);
%     uyVal(i) = uyFun([X(i);Y(i)]);
% end
% 
% figure;
% for i = 1:nX-1
%     for j = 1:nY-1
%         vertX = [X(i,j)+uxVal(i,j); X(i+1,j)+uxVal(i+1,j); X(i+1,j+1)+uxVal(i+1,j+1); X(i,j+1)+uxVal(i,j+1)];
%         vertY = [Y(i,j)+uyVal(i,j); Y(i+1,j)+uyVal(i+1,j); Y(i+1,j+1)+uyVal(i+1,j+1); Y(i,j+1)+uyVal(i,j+1)];
%         patch('Faces',1:4,'Vertices',[vertX,vertY],'EdgeColor','black','FaceColor','none');
%     end
% end
% title('Real solution: u');

%%
fprintf('h: %f e(u)_h0: %f e(u)_h1: %f e(z)_h0: %f e(z)_div: %f e(p)_h0: %f e(Php)_h0: %f\n',msh.h,eu_h0,eu_h1,ez_h0,ez_div,ep_h0,ePhp_h0);
Err = table(eu_h0,eu_h1,ez_h0,ez_div,ep_h0,ePhp_h0);
end