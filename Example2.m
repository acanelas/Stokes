% Analysis of Example 2

a1 = figure;
clear EST;
EST.YC = 0.0;
Example2Zaux;

a2 = figure;
plot([EST.MatNod1(:,1);EST.MatNod1(1,1)],[EST.MatNod1(:,2);EST.MatNod1(1,2)],'k','linewidth',1);
figurita2(EST);

for YC = 0.04:0.04:0.16
    figure(a1);
    EST.YC = YC;
    Example2Zaux;
    
    figure(a2);
    figurita2(EST)
end

figure(a1);
clear EST;
EST.YC = 0.0;
Example2Zaux;

for YC = -0.04:-0.04:-0.24
    figure(a1);
    EST.YC = YC;
    Example2Zaux;

    figure(a2);
    figurita2(EST)
end
