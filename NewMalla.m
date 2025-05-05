function EST = NewMalla(EST)

Lmax = EST.Delta;
Beta = EST.Beta;

NumNod = EST.NumNod2;
NumEle = NumNod;
MatNod = [(0:NumNod-1)',(1:NumNod)'];
MatNod(1,1) = NumNod;
MatEle = [(1:NumNod)',(2:NumNod+1)'];
MatEle(NumNod,2) = 1;
XO = EST.MatNod2;

% Calculo curvaturas de nodos ---------------------------------------------
KNod = zeros(NumNod,1);

for Nod = 1:NumNod
    Ele1 = MatNod(Nod,1);
    Nod1 = MatEle(Ele1,1);
    Nod2 = MatEle(Ele1,2);
    
    X1 = XO(Nod1,:);
    X2 = XO(Nod2,:);
    
    Ele2 = MatNod(Nod,2);
    Nod1 = MatEle(Ele2,1);
    Nod2 = MatEle(Ele2,2);
    
    X3 = XO(Nod1,:);
    X4 = XO(Nod2,:); % X22 = X2;
    
    L1 = norm(X2-X1);
    L2 = norm(X4-X3);
    L3 = norm(X4-X1);
    KNod(Nod) = sqrt((L1+L2+L3)*(L1+L2-L3)*(L3+L1-L2)*(L2+L3-L1))/(L1*L2*L3);
end

% Calculo largos maximos en elementos -------------------------------------
LM = zeros(NumEle,4);

for Ele = 1:NumEle
    Nod1 = MatEle(Ele,1);
    Nod2 = MatEle(Ele,2);
    
    X1 = XO(Nod1,:);
    X2 = XO(Nod2,:);
    
    K1 = KNod(Nod1);
    K2 = KNod(Nod2);
    
    L = Lmax;
    if not(K1==0), L = min(L,Beta/K1); end
    if not(K2==0), L = min(L,Beta/K2); end
    
    LM(Ele,1) = norm(X2-X1); % Longitud del elemento
    LM(Ele,2) = L;           % Longitud maxima en el elemento
    LM(Ele,3:4) = (X2-X1)/LM(Ele,1);
end

% Chequeo necesidad de remallado ------------------------------------------
if any(LM(:,1) > (3/2)*LM(:,2)) || any(LM(:,1) < (2/3)*LM(:,2))
    
% Primera tentativa de malla ----------------------------------------------
XP = XO(1,:);
Ele = 1;
sigo = 1;
Lm = inf; % Largo maximo actual

NumNodo = 1;
MatNodo = zeros(10*NumNod,4);
MatNodo(1,1:3) = [1, XP]; % Elemento y punto

while sigo
    Nod1 = MatEle(Ele,1);
    
    X1 = XO(Nod1,:);
    Lm = min(Lm,LM(Ele,2));
    
    T = LM(Ele,3:4);
    
    a = (X1-XP)*T';        % a puede tener cualquier signo
    h2 = norm(X1-XP)^2 - a^2; % el resultado es positivo
    alpha = sqrt(abs(Lm^2 - h2)) - a;
    
    if (alpha < 0) % puede suceder si se modifico Lm
        alpha = 0;
    end
    if (alpha < LM(Ele,1))
        % Adiciono nodo
        XPn = X1 + alpha*T;
        NumNodo = NumNodo+1;
        MatNodo(NumNodo,:) = [Ele, XPn, norm(XPn-XP)];
        sigo = true;
        XP = XPn;
        Lm = inf;
    else
        % Continuo en el siguiente elemento
        Ele = Ele+1;
        if (Ele > NumEle)
            sigo = false;
            X2 = XO(1,:);
            DL = norm(X2-XP);
            Long = sum(MatNodo(:,4));
            gamma = (Long+DL)/(Long+Lm);
            MatNodo = MatNodo(1:NumNodo,:);
            MatNodo(:,4) = gamma*MatNodo(:,4);
        end
    end
end

% Segunda tentativa de malla ----------------------------------------------
XP = XO(1,:);
Ele = 1;
sigo = 1;

NumNodoAnterior = NumNodo;
NumNodo = 1;
MatNodo(1,1:3) = [1, XP]; % Elemento y punto

while sigo
    Nod1 = MatEle(Ele,1);
    
    X1 = XO(Nod1,:);
    Lm = MatNodo(NumNodo+1,4);
    
    T = LM(Ele,3:4);
    
    a = (X1-XP)*T';        % a puede tener cualquier signo
    h2 = norm(X1-XP)^2 - a^2; % el resultado es positivo
    alpha = sqrt(abs(Lm^2 - h2)) - a;
    
    if (alpha < 0) % puede suceder si se modifico Lm
        alpha = 0;
    end
    if (alpha < LM(Ele,1))
        % Adiciono nodo
        XPn = X1 + alpha*T;
        NumNodo = NumNodo+1;
        MatNodo(NumNodo,1:3) = [Ele, XPn];
        sigo = true;
        XP = XPn;
    else
        % Continuo en el siguiente elemento
        Ele = Ele+1;
    end
    if (NumNodo >= NumNodoAnterior)
        sigo = false;
    end
end

% Tercera tentativa de malla (Suavizo) ------------------------------------
for Nod = 1:NumNodo
    Ele = MatNodo(Nod,1);
    EleA = Ele-1;
    if (EleA<1), EleA = NumNod; end
    EleB = Ele+1;
    if (EleB>NumNod), EleB = 1; end
    
    X  = MatNodo(Nod,2:3);
    
    X1 = XO(MatEle(EleA,1),:); S1 = 0;
    X2 = XO(MatEle(EleA,2),:); S2 = LM(EleA,1);
    X3 = XO(MatEle(EleB,1),:); S3 = S2 + LM(Ele,1); 
    X4 = XO(MatEle(EleB,2),:); S4 = S3 + LM(EleB,1);
    
    S = S2 + norm(X-X2);
    
    N1 = (S-S2)*(S-S3)*(S-S4)/(S1-S2)/(S1-S3)/(S1-S4);
    N2 = (S-S1)*(S-S3)*(S-S4)/(S2-S1)/(S2-S3)/(S2-S4);
    N3 = (S-S1)*(S-S2)*(S-S4)/(S3-S1)/(S3-S2)/(S3-S4);
    N4 = (S-S1)*(S-S2)*(S-S3)/(S4-S1)/(S4-S2)/(S4-S3);
    
    XN = X1*N1+X2*N2+X3*N3+X4*N4;
    MatNodo(Nod,1:3) = [Ele, XN];
end

% interpolacion de variables ----------------------------------------------
XU = zeros(2*NumNodo,1);
XLAM = EST.x0(2*EST.NumNod1+3*EST.NumNod2+1);

for Nod = 1:NumNodo
    Ele = MatNodo(Nod,1);
    XP = MatNodo(Nod,2:3);
    
    Nod1 = MatEle(Ele,1);
    Nod2 = MatEle(Ele,2);
    
    X1 = XO(Nod1,:);
    X2 = XO(Nod2,:);
    
    LE = norm(X2-X1);
    N1 = 1 - norm(XP-X1)/LE;
    N2 = 1 - norm(XP-X2)/LE;
    
    UX1 = EST.XU(2*Nod1-1);
    UY1 = EST.XU(2*Nod1);
    UX2 = EST.XU(2*Nod2-1);
    UY2 = EST.XU(2*Nod2);
    
    XU(2*Nod-1) = N1*UX1 + N2*UX2;
    XU(2*Nod)   = N1*UY1 + N2*UY2;
end

EST.NumNod2 = NumNodo;
EST.NumEle2 = NumNodo;
EST.MatNod2 = MatNodo(:,2:3);
EST.MatEle2 = [(1:NumNodo)', (2:NumNodo+1)'];
EST.MatEle2(NumNodo,2) = 1;

EST.XU = XU;

% Velocidades -------------------------------------------------------------
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

% x0 ----------------------------------------------------------------------
NE1 = EST.NumEle1; NE2 = EST.NumEle2;
NumVar = 2*NE1 + 3*NE2 + 1;

EST.x0 = zeros(2*EST.NumNod1+3*EST.NumNod2+1,1);
EST.x0(1:2*NE1)             = EST.XP;
EST.x0(2*NE1+1:2*NE1+2*NE2) = EST.XU;
EST.x0(NumVar) = XLAM;

end
