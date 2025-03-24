function myplot(X,Y,xlbl,ylbl)


figr = figure;
axfigr = axes('Parent',figr);
hold(axfigr,'on');


plot(X,Y,'MarkerSize',10,'Marker','.','LineWidth',1);

box(axfigr,'on');
axis(axfigr,'square');
hold(axfigr,'off');
set(axfigr,'FontSize',15,'LineWidth',1,'XMinorTick','on','YMinorTick','on');
ylabel(ylbl,'FontSize',30,'Interpreter','latex','Rotation',0);
xlabel(xlbl,'FontSize',30,'Interpreter','latex');

end