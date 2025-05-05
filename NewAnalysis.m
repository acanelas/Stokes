function EST = NewAnalysis(EST)

% Velocidades -------------------------------------------------------------
EST.NumNod1 = size(EST.MatNod1,1);
EST.NumNod2 = size(EST.MatNod2,1);
EST.NumEle1 = size(EST.MatNod1,1);
EST.NumEle2 = size(EST.MatNod2,1);

EST.VZ = zeros(EST.NumNod2,2);
for Ele = 1:EST.NumEle2
    Nod1 = EST.MatEle2(Ele,1);
    Nod2 = EST.MatEle2(Ele,2);
    X1 = EST.MatNod2(Nod1,1);
    Y1 = EST.MatNod2(Nod1,2);
    X2 = EST.MatNod2(Nod2,1);
    Y2 = EST.MatNod2(Nod2,2);
    V = [(Y2-Y1), -(X2-X1)];
    V = V/norm(V);
    EST.VZ(Nod1,:) = EST.VZ(Nod1,:) + V;
    EST.VZ(Nod2,:) = EST.VZ(Nod2,:) + V;
end
for Nod = 1:EST.NumNod2
    EST.VZ(Nod,:) = EST.VZ(Nod,:)/norm(EST.VZ(Nod,:));
end

if EST.REM
    EST = NewMalla(EST);
end

% Calculo presion interior normal -----------------------------------------
EST.PN = zeros(2*EST.NumNod1,1);
Normal = zeros(EST.NumNod1,2);
for Ele = 1:EST.NumEle1
    Nod1 = EST.MatEle1(Ele,1);
    Nod2 = EST.MatEle1(Ele,2);
    X1 = EST.MatNod1(Nod1,1);
    Y1 = EST.MatNod1(Nod1,2);
    X2 = EST.MatNod1(Nod2,1);
    Y2 = EST.MatNod1(Nod2,2);
    V = [(Y2-Y1), -(X2-X1)];
    V = V/norm(V);
    Normal(Nod1,:) = Normal(Nod1,:) + V;
    Normal(Nod2,:) = Normal(Nod2,:) + V;
end
for Nod = 1:EST.NumNod1
    Normal(Nod,:) = Normal(Nod,:)/norm(Normal(Nod,:));
    EST.PN(2*Nod-1) = Normal(Nod,1);
    EST.PN(2*Nod)   = Normal(Nod,2);
end
    
% Iteracion principal -----------------------------------------------------
F0 = EST.fun(EST);
NF0 = norm(F0);

fprintf('----------------------------------------------------\n');
fprintf('  ITER | ITERI | NumNod |    PASO    |      NF0     \n');
fprintf('%6i |%6i |%7i | %10.6f | %12.4e\n',0,0,EST.NumNod2,0,NF0);

ITER = 1;
while ((NF0>EST.TOL)&&(ITER<EST.MaxIt)||(ITER<=EST.MinIt))
    gradF0 = EST.gfun(EST);
    
    DX = gradF0\(-F0);
    B = max(abs(DX(2*EST.NumNod1+2*EST.NumNod2+1:2*EST.NumNod1+3*EST.NumNod2)));
    A = min(B,EST.PasoMax);
    AB = A/B;
    DX = DX*AB;
    
    x0 = EST.x0;
    x1 = x0 + DX;
    EST.x0 = x1;
    F1 = EST.fun(EST);
    NF1 = norm(F1);
    
    ITERI = 0;
    while ((NF1 > NF0)&&(ITERI<10))
        ITERI = ITERI + 1;
        AB = 0.5*AB;
        DX = DX*0.5;
        x1 = x0 + DX;
        EST.x0 = x1;
        F1 = EST.fun(EST);
        NF1 = norm(F1);
    end
    
    % Actualizo datos
    EST = NewDatos(EST);
    
    if EST.REM
        EST = NewMalla(EST);
    end

    F0 = EST.fun(EST);
    NF0 = norm(F0);

    fprintf('%6i |%6i |%7i | %10.6f | %12.4e\n',ITER,ITERI,EST.NumNod2,AB,NF0);

    figurita(EST);

    ITER = ITER+1;
end

disp('Norm of F:');
disp(NF0);
