% Analysis of Example 4

a1 = figure;
clear EST;
EST.YC = 0.0;
Example4Zaux;

a2 = figure;
plot([EST.MatNod1(:,1);EST.MatNod1(1,1)],[EST.MatNod1(:,2);EST.MatNod1(1,2)],'k','linewidth',1);
figurita2(EST);

for YC = -0.002:-0.002:-0.04
    figure(a1);
    EST.YC = YC;
    Example2Zaux;
    
    % figure(a2);
    % figurita2(EST)
end

for YC = -0.04:-0.04:-0.40
    figure(a1);
    EST.YC = YC;
    Example4Zaux;

    figure(a2);
    figurita2(EST)
end
