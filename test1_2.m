%% Parameter Setting
domn = [0,1,0,1]; nSub = 100;
EgDire = [2,3,4]; EgNeum = 1;
lmd = 1; mu = 1; alph = 1; c0 = 0; k = 1;
vu = 0.1; vps = 10.^(-1:6);
test = 'test1'; method = 'BR'; index = '1';

ux = sprintf('%.1g*sin(pi*x)*(1-cos(pi*y))',vu);
uy = sprintf('%.1g*(1-cos(2*pi*x))*sin(pi/2*y)',vu);
p = 'vp*cos(pi*x)*cos(pi*y)';

%% Main 
Domn = struct('domn',domn,'nSub',nSub,'EgDire',EgDire,'EgNeum',EgNeum);
Parm = struct('lmd',lmd,'mu',mu,'alph',alph,'c0',c0,'k',k);
Sol = struct('ux',ux,'uy',uy,'p',[]);

eEvl = table('Size',[length(vps),7],...
    'VariableTypes',repmat("double",[1,7]),...
    'VariableNames',{'vp','eu_h0','eu_h1','ez_h0','ez_div','ep_h0','ePhp_h0'});

for i = 1:length(vps)
    vp = vps(i); Sol.p = strrep(p,'vp',num2str(vp));
%     Err = Biot_LR(Domn,Parm,Sol);
    eval(['Err = Biot_', method, '(Domn,Parm,Sol,[]);'])
    eEvl(i,1) = {vp}; eEvl(i,2:7) = Err;
end

eEvl.Properties.Description = sprintf(['======== %s ========\n ',...
    'ux = %s\n uy= %s\n p = %s\n method: %s index: %s \n ',...
    'nSub = %u lmd = %.1g mu = %.1g alph = %.1g c0 = %.1g k = %.1g vu = %.1g'],...
    test,ux,uy,p,method,index,nSub,lmd,mu,alph,c0,k,vu);

disp(eEvl.Properties.Description);
disp(eEvl);

%% Output
file = [test, '_result']; 
name = ['pRob_',method,'_',index];
eval([name,' = eEvl;']);
save([file,'\',name,'.mat'],name);