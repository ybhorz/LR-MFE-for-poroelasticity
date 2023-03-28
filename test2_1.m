%% Parameter Setting
domn = [0,1,0,1]; nSubs = [8,16,32,64,128];
EgDire = [2,3,4]; EgNeum = 1;
lmd = 1e6; mu = 1; alph = 1; c0 = 0; k = 1;
test = 'test2'; method = 'P1'; index = '1';

ux = '0.1*(sin(2*pi*y)*(-1+cos(2*pi*x))+1/lmd*sin(pi*x)*sin(pi*y))';
uy = '0.1*(sin(2*pi*x)*(1-cos(2*pi*y))+1/lmd*sin(pi*x)*sin(pi*y))';
p = 'cos(pi*x)*(1-cos(pi*y))';

%% Main 
Domn = struct('domn',domn,'nSub',[],'EgDire',EgDire,'EgNeum',EgNeum);
Parm = struct('lmd',lmd,'mu',mu,'alph',alph,'c0',c0,'k',k);
Sol = struct('ux',ux,'uy',uy,'p',p);

eRate = table('Size',[length(nSubs),13],...
    'VariableTypes',repmat("double",[1,13]),...
    'VariableNames',{'nSub','eu_h0','ru_h0','eu_h1','ru_h1','ez_h0','rz_h0',...
    'ez_div','rz_div','ep_h0','rp_h0','ePhp_h0','rPhp_h0'});

for i = 1:length(nSubs)
    nSub = nSubs(i); Domn.nSub = nSub;
%     Err = Biot_LR(Domn,Parm,Sol);
    eval(['Err = Biot_', method, '(Domn,Parm,Sol,[]);'])
    eRate(i,1) = {nSub}; eRate(i,2:2:13) = Err;
    if i >= 2
        for j = 3:2:13
            eRate{i,j} = log(eRate{i-1,j-1}/eRate{i,j-1})/log(2);
        end
    end
end

eRate.Properties.Description = sprintf(['======== %s ========\n ',...
    'ux = %s\n uy= %s\n p = %s\n method: %s index: %s \n ',...
    'lmd = %.1g mu = %.1g alph = %.1g c0 = %.1g k = %.1g'],...
    test,ux,uy,p,method,index,lmd,mu,alph,c0,k);

disp(eRate.Properties.Description);
disp(eRate);
%% Output
file = [test, '_result']; 
name = ['eRate_',method,'_',index];
eval([name,' = eRate;']);
save([file,'\',name,'.mat'],name);