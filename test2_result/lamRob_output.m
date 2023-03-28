load('lmdRob_LR_1.mat');
figure;
semilogx(lmdRob_LR_1,"lmd","eu_h1","LineWidth",1,"Marker","o"); hold on;
semilogx(lmdRob_LR_1,"lmd","ez_h0","LineWidth",1,"Marker","square");
semilogx(lmdRob_LR_1,"lmd","ez_div","LineWidth",1,"Marker","diamond");
semilogx(lmdRob_LR_1,"lmd","ep_h0","LineWidth",1,"Marker","^"); hold off;
axis([5e-1,2e7,-0.05,1.2]);
xticks(10.^(0:7)); yticks(0:0.2:1.2);
xlabel('$\lambda$','Interpreter','latex');
ylabel('error');
legend('$|||\mathbf{u} - \mathbf{u}_h|||$', '$||\mathbf{z}-\mathbf{z}_h||$',...
    '$||\nabla\cdot (\mathbf{z} - \mathbf{z}_h)||$', '$||p-p_h||$',...
    'Interpreter','latex','Location','northwest');
title('Error evolution of LR-MFE method');
%%
load('lmdRob_P1_1.mat');
figure;
semilogx(lmdRob_P1_1,"lmd","eu_h1","LineWidth",1,"Marker","o"); hold on;
semilogx(lmdRob_P1_1,"lmd","ez_h0","LineWidth",1,"Marker","square");
semilogx(lmdRob_P1_1,"lmd","ez_div","LineWidth",1,"Marker","diamond");
semilogx(lmdRob_P1_1,"lmd","ep_h0","LineWidth",1,"Marker","^"); hold off;
axis([5e-1,2e7,-0.05,1.2]);
xticks(10.^(0:7)); yticks(0:0.2:1.2);
xlabel('$\lambda$','Interpreter','latex');
ylabel('error');
legend('$|||\mathbf{u} - \mathbf{u}_h|||$', '$||\mathbf{z}-\mathbf{z}_h||$',...
    '$||\nabla\cdot (\mathbf{z} - \mathbf{z}_h)||$', '$||p-p_h||$',...
    'Interpreter','latex','Location','northwest');
title('Error evolution of CG-MFE method');