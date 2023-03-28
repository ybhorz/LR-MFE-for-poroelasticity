figure;
rectangle('Position',[0 0 1 1]);
axis([-0.25,1.25,-0.25,1.25]); axis off;
text(0.5,1.1,{'{\boldmath$\tilde{\sigma}$} $\mathbf{n} = \mathbf{t}_N$', '$\mathbf{z}\cdot \mathbf{n} = 0$'},'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','bottom');
text(1.1,0.5,{'{\boldmath$\tilde{\sigma}$} $\mathbf{n} = \mathbf{0}$', '$\mathbf{z}\cdot \mathbf{n} = 0$'},'Interpreter','latex','HorizontalAlignment','left','VerticalAlignment','middle');
text(0.5,-0.1,{'{\boldmath$\tilde{\sigma}$} $\mathbf{n} = \mathbf{0}$', '$\mathbf{z}\cdot \mathbf{n} = 0$'},'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','top');
text(-0.1,0.5,{'$\mathbf{u} = \mathbf{0}$', '$\mathbf{z}\cdot \mathbf{n} = 0$'},'Interpreter','latex','HorizontalAlignment','right','VerticalAlignment','middle');
text(0.5,0.5,{'$\Omega$'},'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','middle');
