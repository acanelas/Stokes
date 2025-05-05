function EST = NewDatos(EST)

% Actualizo ---------------------------------------------------------------
x1 = EST.x0(2*EST.NumNod1+2*EST.NumNod2+1:2*EST.NumNod1+3*EST.NumNod2);
EST.x0(2*EST.NumNod1+2*EST.NumNod2+1:2*EST.NumNod1+3*EST.NumNod2) = 0;

for Nod = 1:EST.NumNod2
    EST.MatNod2(Nod,1) = EST.MatNod2(Nod,1) + x1(Nod)*EST.VZ(Nod,1);
    EST.MatNod2(Nod,2) = EST.MatNod2(Nod,2) + x1(Nod)*EST.VZ(Nod,2);
end

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

% Otros -------------------------------------------------------------------
EST.XP = EST.x0(1:2*EST.NumNod1);
EST.XU = EST.x0(2*EST.NumNod1+1:2*EST.NumNod1+2*EST.NumNod2);
