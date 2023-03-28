%% tri-subdivision: mesh2D_h0 (h), mesh2D_h1 (h/2), mesh2D_h2 (h/4)
clear; close all;
% domn = [0,1,0,1];
% load("mesh2D_[0,1,0,1]\mesh2D_h1.mat");
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
%% unifrom-triangulation
domn = [0,1,0,1];
nSub = 16;
[node,elem,edge,conn,elemCn] = genTri(domn,nSub);
nNode = size(node,2); nElem = size(elem,2); nEdge = edge.nConn; nConn = conn.nConn;
%% Real solution
lmd = 1; mu = 1; alph = 1; c0 = 0; k = 1e-8;

syms x y;
var = [x;y];
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

% Neuman boundary condition.
Fv = SLF(Fcn(var,-1e-1),1,d0uy,"edge",find(msh.edge.type==1));
%% Assmeble matrix
Stiff11 = assemble(msh,Uh,Uh,Auv,SLF.empty,"Integral",GInt('D2P1')); Stiff11 = sparse(Stiff11);
Stiff12 = sparse(Uh.nNode,Zh.nNode);
Stiff13 = assemble(msh,Ph,Uh,Bpv,SLF.empty,"Integral",GInt('D2P1')); Stiff13 = sparse(Stiff13);
[~,Load1] = assemble(msh,FE.empty,Uh,DLF.empty,Fv,"Integral",GInt('D2P4'));

Stiff21 = sparse(Zh.nNode,Uh.nNode);
Stiff22 = assemble(msh,Zh,Zh,Azw,SLF.empty,"MatrixType","sparse","Integral",GInt('D2P2')); Stiff22 = sparse(Stiff22);
Stiff23 = assemble(msh,Ph,Zh,Bpw,SLF.empty,"Integral",GInt('D2P1')); Stiff23 = sparse(Stiff23);
Load2 = zeros([Zh.nNode,1]);

Stiff31 = assemble(msh,Uh,Ph,Buq,SLF.empty,"Integral",GInt('D2P1')); Stiff31 = sparse(Stiff31);
Stiff32 = assemble(msh,Zh,Ph,Bzq,SLF.empty,"Integral",GInt('D2P1')); Stiff32 = sparse(Stiff32);
Stiff33 = sparse(Ph.nNode,Ph.nNode);
Load3 = zeros([Ph.nNode,1]);
%% Dirichlet boundary condition
% type 1: y = 1
% type 2: x = 1
% type 3: y = 0 
% type 4: x = 0

for iEdge = 1:Uh.nEdge 
    for iNode = 1:2
        I = Uh.edge.node(iNode,iEdge);
        coord = Uh.node(:,I);
        if ismember(Uh.edge.type(iEdge),4)
            Stiff11(I,:) = 0;
            Stiff11(I,I) = 1;
            Stiff13(I,:) = 0;
            Load1(I) = 0;
        end
    end
    for iNode = 3:4
        I = Uh.edge.node(iNode,iEdge);
        coord = Uh.node(:,I);
        if ismember(Uh.edge.type(iEdge),4)
            Stiff11(I,:) = 0;
            Stiff11(I,I) = 1;
            Stiff13(I,:) = 0;
            Load1(I) = 0;
        end
    end
end

for iEdge = 1:Zh.nEdge
    I = Zh.edge.node(iEdge);
    Stiff22(I,:) = 0;
    Stiff22(I,I) = 1;
    Stiff23(I,:) = 0;
    Load2(I) = 0;
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

uPlot(domn,msh,Uh,uh);
pPlot(msh,Ph,ph);
%% Plot: p
function pPlot(msh,Ph,ph)
    figure; view(3);
    for iElem = 1:msh.nElem
        phVal = repmat(ph(Ph.elem(:,iElem)),[3,1]);
        patch('Faces',1:3,'Vertices',[msh.node(:,msh.elem(:,iElem))',phVal],'FaceVertexCData',phVal,'FaceColor','interp');
    end
    title('Pressure');
    view(-140,30)
end
%% Plot: u
function uPlot(domn,msh,Uh,uh)
    nUhxNodes = 1:msh.nNode;
    nUhyNodes = msh.nNode + (1:msh.nNode);
    
    [X,Y] = meshgrid(domn(1):msh.h:domn(2),domn(3):msh.h:domn(4));
    nX = size(X,2);
    nY = size(Y,1);
    
    uhx = uh(nUhxNodes);
    uhy = uh(nUhyNodes);
    uhxNode = Uh.node(:,nUhxNodes);
    uhyNode = Uh.node(:,nUhyNodes);
    uhxVal = griddata(uhxNode(1,:)',uhxNode(2,:)',uhx,X,Y,'cubic');
    uhyVal = griddata(uhyNode(1,:)',uhyNode(2,:)',uhy,X,Y,'cubic');
    
    figure;
    for i = 1:nX-1
        for j = 1:nY-1
            vertX = [X(i,j)+uhxVal(i,j); X(i+1,j)+uhxVal(i+1,j); X(i+1,j+1)+uhxVal(i+1,j+1); X(i,j+1)+uhxVal(i,j+1)];
            vertY = [Y(i,j)+uhyVal(i,j); Y(i+1,j)+uhyVal(i+1,j); Y(i+1,j+1)+uhyVal(i+1,j+1); Y(i,j+1)+uhyVal(i,j+1)];
            patch('Faces',1:4,'Vertices',[vertX,vertY],'EdgeColor','black','FaceColor','none');
        end
    end
    axis([0,1.1,-0.1,1])
    xticks(0:0.2:1); yticks(0:0.2:1);
    title('Deformation');
end