%% Parameter Setting
domn = [0,1,0,1]; nSub = 32;
EgDire = [2,3,4]; EgNeum = 1;
lmds = 10.^(0:7); 
mu = 1; alph = 1; c0 = 0; k = 1;
test = 'test2'; method = 'LR'; index = '1';

ux = '0.1*(sin(2*pi*y)*(-1+cos(2*pi*x))+1/lmd*sin(pi*x)*sin(pi*y))';
uy = '0.1*(sin(2*pi*x)*(1-cos(2*pi*y))+1/lmd*sin(pi*x)*sin(pi*y))';
p = 'cos(pi*x)*(1-cos(pi*y))';

%% Main 
Domn = struct('domn',domn,'nSub',nSub,'EgDire',EgDire,'EgNeum',EgNeum);
Parm = struct('lmd',[],'mu',mu,'alph',alph,'c0',c0,'k',k);
Sol = struct('ux',[],'uy',[],'p',p);

eEvl = table('Size',[length(lmds),7],...
    'VariableTypes',repmat("double",[1,7]),...
    'VariableNames',{'lmd','eu_h0','eu_h1','ez_h0','ez_div','ep_h0','ePhp_h0'});

for i = 1:length(lmds)
    lmd = lmds(i); Parm.lmd = lmd;
    Sol.ux = strrep(ux,'lmd',num2str(lmd)); 
    Sol.uy = strrep(uy,'lmd',num2str(lmd));
%     Err = Biot_P1(Domn,Parm,Sol); % Parameter
    eval(['Err = Biot_',method,'(Domn,Parm,Sol,[]);']);
    eEvl(i,1) = {lmd}; eEvl(i,2:end) = Err;
end

%%
eEvl.Properties.Description = sprintf(['======== %s ========\n ',...
    'ux = %s\n uy= %s\n p = %s\n method: %s index: %s \n ',...
    'nSub = %u mu = %.1g alph = %.1g c0 = %.1g k = %.1g'],...
    test,ux,uy,p,method,index,nSub,mu,alph,c0,k);

disp(eEvl.Properties.Description);
disp(eEvl);

%% Output
file = [test, '_result']; 
name = ['lmdRob_',method,'_',index];
eval([name,' = eEvl;']);
save([file,'\',name,'.mat'],name);