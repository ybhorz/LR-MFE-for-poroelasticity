%% tri-subdivision: mesh2D_h0 (h), mesh2D_h1 (h/2), mesh2D_h2 (h/4)
clear; close all;
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

% Neuman boundary condition.
Fv1 = SLF(Fcn(var,-1e-1),1,d0uy,"edge",find(msh.edge.type==1));
FvR = SLF(Fcn(var,-1e-1),1,d0uy,"edge",find(msh.edge.type==1));
%% Assmeble matrix
Stiff11 = assemble(msh,Uh1,Uh1,Auv1,SLF.empty,"Integral",GInt('D2P1')); Stiff11 = sparse(Stiff11);
Stiff12 = assemble(msh,UhR,Uh1,AuRv1,SLF.empty,"MatrixType","sparse","Integral",GInt('D2P1'));
Stiff13 = sparse(Uh1.nNode,Zh.nNode);
Stiff14 = assemble(msh,Ph,Uh1,Bpv1,SLF.empty,"Integral",GInt('D2P1')); Stiff14 = sparse(Stiff14);
[~,Load1] = assemble(msh,FE.empty,Uh1,DLF.empty,Fv1,"Integral",GInt('D2P4'));

Stiff21 = assemble(msh,Uh1,UhR,Au1vR,SLF.empty,"MatrixType","sparse","Integral",GInt('D2P1'));
Stiff22 = assemble(msh,UhR,UhR,AuvR,SLF.empty,"MatrixType","sparse","Integral",GInt('D2P2'));
Stiff23 = sparse(UhR.nNode,Zh.nNode);
Stiff24 = assemble(msh,Ph,UhR,BpvR,SLF.empty,"Integral",GInt('D2P1')); Stiff24 = sparse(Stiff24);
[~,Load2] = assemble(msh,FE.empty,UhR,DLF.empty,FvR,"Integral",GInt('D2P4'));

Stiff31 = sparse(Zh.nNode,Uh1.nNode);
Stiff32 = sparse(Zh.nNode,UhR.nNode);
Stiff33 = assemble(msh,Zh,Zh,Azw,SLF.empty,"MatrixType","sparse","Integral",GInt('D2P2')); Stiff33 = sparse(Stiff33);
Stiff34 = assemble(msh,Ph,Zh,Bpw,SLF.empty,"Integral",GInt('D2P1')); Stiff34 = sparse(Stiff34);
Load3 = zeros([Zh.nNode,1]);

Stiff41 = assemble(msh,Uh1,Ph,Buq1,SLF.empty,"Integral",GInt('D2P1')); Stiff41 = sparse(Stiff41);
Stiff42 = assemble(msh,UhR,Ph,BuqR,SLF.empty,"Integral",GInt('D2P1')); Stiff42 = sparse(Stiff42);
Stiff43 = assemble(msh,Zh,Ph,Bzq,SLF.empty,"Integral",GInt('D2P1')); Stiff43 = sparse(Stiff43);
Stiff44 = sparse(Ph.nNode,Ph.nNode);
Load4 = zeros([Ph.nNode,1]);
%% Dirichlet boundary condition
% type 1: y = 1
% type 2: x = 1
% type 3: y = 0 
% type 4: x = 0

for iEdge = 1:Uh1.nEdge 
    for iNode = 1:2
        I = Uh1.edge.node(iNode,iEdge);
        coord = Uh1.node(:,I);
        if ismember(Uh1.edge.type(iEdge),4)
            Stiff11(I,:) = 0;
            Stiff11(I,I) = 1;
            Stiff12(I,:) = 0;
            Stiff14(I,:) = 0;
            Load1(I) = 0;
        end
    end
    for iNode = 3:4
        I = Uh1.edge.node(iNode,iEdge);
        coord = Uh1.node(:,I);
        if ismember(Uh1.edge.type(iEdge),4)
            Stiff11(I,:) = 0;
            Stiff11(I,I) = 1;
            Stiff12(I,:) = 0;
            Stiff14(I,:) = 0;
            Load1(I) = 0;
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

for iEdge = 1:Zh.nEdge
    I = Zh.edge.node(iEdge);
    Stiff33(I,:) = 0;
    Stiff33(I,I) = 1;
    Stiff34(I,:) = 0;
    Load3(I) = 0;
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

uPlot(domn,msh,Uh1,uh1);
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
    axis([0,1.1,-0.15,1])
    xticks(0:0.2:1); yticks(0:0.2:1);
    title('Deformation');
end