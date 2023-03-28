figure;
rectangle('Position',[0 0 1 1]);
axis([-0.75,1.75,-0.2,1.2]); axis off;
text(0.5,1.05,'$u_x = \frac{\partial u_y}{\partial y} = p =0$','Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','bottom');
text(1.05,0.5,'$\frac{\partial u_x}{\partial x} = u_y = p =0$','Interpreter','latex','HorizontalAlignment','left','VerticalAlignment','middle');
text(0.5,-0.05,'$u_x = \frac{\partial u_y}{\partial y} = p =0$','Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','top');
text(-0.05,0.5,'$\frac{\partial u_x}{\partial x} = u_y = p =0$','Interpreter','latex','HorizontalAlignment','right','VerticalAlignment','middle');
text(0.26,0.26,{'$(x_0,y_0)$'},'Interpreter','latex','HorizontalAlignment','left','VerticalAlignment','bottom');
text(0.25,0.25,{'$\bullet$'},'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','middle');