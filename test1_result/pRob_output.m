load('pRob_LR_1.mat'); pRob_LR_1 = pRob_LR_1(2:end,:);
figure;
loglog(pRob_LR_1,"vp","eu_h1","LineWidth",1,"Marker","o"); hold on;
loglog(pRob_LR_1,"vp","ez_h0","LineWidth",1,"Marker","square");
loglog(pRob_LR_1,"vp","ez_div","LineWidth",1,"Marker","diamond");
loglog(pRob_LR_1,"vp","ep_h0","LineWidth",1,"Marker","^"); hold off;
axis([5e-1,2e6,2e-3,5e5]);
xticks(10.^(0:6)); yticks(10.^(-2:5));
xlabel('$\nu_p$','Interpreter','latex');
ylabel('error');
legend('$|||\mathbf{u} - \mathbf{u}_h|||$', '$||\mathbf{z}-\mathbf{z}_h||$',...
    '$||\nabla\cdot (\mathbf{z} - \mathbf{z}_h)||$', '$||p-p_h||$',...
    'Interpreter','latex','Location','northwest');
title('Error evolution of LR-MFE method');
%%
load('pRob_BR_1.mat'); pRob_BR_1 = pRob_BR_1(2:end,:);
figure;
loglog(pRob_BR_1,"vp","eu_h1","LineWidth",1,"Marker","o"); hold on;
loglog(pRob_BR_1,"vp","ez_h0","LineWidth",1,"Marker","square");
loglog(pRob_BR_1,"vp","ez_div","LineWidth",1,"Marker","diamond");
loglog(pRob_BR_1,"vp","ep_h0","LineWidth",1,"Marker","^"); hold off;
axis([5e-1,2e6,2e-3,5e5]);
xticks(10.^(0:6)); yticks(10.^(-2:5));
xlabel('$\nu_p$','Interpreter','latex');
ylabel('error');
legend('$|||\mathbf{u} - \mathbf{u}_h|||$', '$||\mathbf{z}-\mathbf{z}_h||$',...
    '$||\nabla\cdot (\mathbf{z} - \mathbf{z}_h)||$', '$||p-p_h||$',...
    'Interpreter','latex','Location','northwest');
title('Error evolution of BR-MFE method');