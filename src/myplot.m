function myplot(X,Y,xlbl,ylbl)


figr = figure('Position',[100,100,512,512]);
axfigr = axes('Parent',figr);
hold(axfigr,'on');


plot(X,Y,'MarkerSize',10,'Marker','.','LineWidth',1);

box(axfigr,'on');
axis(axfigr,'square');
hold(axfigr,'off');
set(axfigr,'FontSize',15,'LineWidth',2,'TickLength',[0.02 0.02],'XMinorTick',...
    'on','YMinorTick','on');
ylabel(ylbl,'FontSize',30,'Interpreter','latex','Rotation',0);
xlabel(xlbl,'FontSize',30,'Interpreter','latex');

lengthX = abs(X(end) - X(1));
Xlowerlim = min(X) - lengthX/32;
Xupperlim = max(X) + lengthX/32;
xlim(axfigr,[Xlowerlim Xupperlim]);

lengthY = abs(Y(end) - Y(1));
Ylowerlim = min(Y) - lengthY/32;
Yupperlim = max(Y) + lengthY/32;
ylim(axfigr,[Ylowerlim Yupperlim]);

end