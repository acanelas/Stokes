function figurita(EST)

plot([EST.MatNod1(:,1);EST.MatNod1(1,1)],[EST.MatNod1(:,2);EST.MatNod1(1,2)],'k','linewidth',1,'Marker','.','MarkerSize',8);
hold on;
plot([EST.MatNod2(:,1);EST.MatNod2(1,1)],[EST.MatNod2(:,2);EST.MatNod2(1,2)],'r','linewidth',1,'Marker','.','MarkerSize',8);

hold off;
axis equal;
pause(0.1)
