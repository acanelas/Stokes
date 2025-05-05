function figurita2(EST)

hold on;
plot([EST.MatNod2(:,1);EST.MatNod2(1,1)],[EST.MatNod2(:,2);EST.MatNod2(1,2)],'r','linewidth',1);

axis equal;
pause(0.1)
