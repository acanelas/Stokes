function G = funFU(EST)

% Unmodified function (F)

EST.NumVar = 2*EST.NumNod1 + 3*EST.NumNod2 + 1; % Nr. of Variables
EST.NumRes = 2*EST.NumNod1 + 3*EST.NumNod2 + 1; % Nr. of Constraints

% Constantes --------------------------------------------------------------

NG = 6;

CD0 = 0.0d0;
CD1 = 1.0d0;
CD2 = 2.0d0;
CD3 = 3.0d0;
CD4 = 4.0d0;
CD5 = 5.0d-1;
CD6 = 1.5d0;

CI3 = CD1/CD3;
CI6 = CD5/CD3;

XX = 1;
YY = 2;

Ge = EST.Mu;
ve = CD5;

MC1 = CD1-ve;
MC2 = CD1-CD2*ve;
MC3 = CD3-CD4*ve;
MC4 = -CD1/(CD4*pi*MC1);
MC5 = MC4/(CD2*Ge);
MC6 = MC3*MC5;
MC7 = MC2*MC4;

TolRa = 1.0d-6;

CM1 = MC4;
CM2 = CD2*CM1;
CM3 = MC5;
CM4 = CD5*CM3;
CM5 = MC6;
CM6 = CD5*CM5;
CM7 = MC7;

% Codigo del programa -----------------------------------------------------

WI = zeros(1,6);
XI = zeros(1,6);

XI(1) = -0.932469514203152d0;	XI(6) = -XI(1);
XI(2) = -0.661209386466265d0;	XI(5) = -XI(2);
XI(3) = -0.238619186083197d0;	XI(4) = -XI(3);
WI(1) =  0.171324492379170d0;	WI(6) =  WI(1);
WI(2) =  0.360761573048139d0;	WI(5) =  WI(2);
WI(3) =  0.467913934572691d0;	WI(4) =  WI(3);

X = EST.x0;

MatXP = 1:2*EST.NumNod1;
MatXU = 2*EST.NumNod1+1:2*EST.NumNod1+2*EST.NumNod2;
MatXT = 2*EST.NumNod1+2*EST.NumNod2+1:2*EST.NumNod1+3*EST.NumNod2;
MatXL = 2*EST.NumNod1+3*EST.NumNod2+1; % Variable lambda adicionada al final

ResLAM = EST.NumRes; % Restriccion de variable lambda

G = zeros(EST.NumRes,1);

for NodF = 1:EST.NumNod1
    
    % Nodo Fuente ---------------------------------------------------------
    
    XF = EST.MatNod1(NodF,XX);
    YF = EST.MatNod1(NodF,YY);
    
    XUFX = EST.XG(2*NodF-1);
    XUFY = EST.XG(2*NodF);
    
    ResF1 = 2*NodF-1;
    ResF2 = 2*NodF;
    
    for Ele = 1:EST.NumEle1
        
        Nod1 = EST.MatEle1(Ele,1);
        Nod2 = EST.MatEle1(Ele,2);
        
        VarP1X = MatXP(2*Nod1-1);
        VarP1Y = MatXP(2*Nod1);
        VarP2X = MatXP(2*Nod2-1);
        VarP2Y = MatXP(2*Nod2);
        
        XP1X = X(VarP1X);
        XP1Y = X(VarP1Y);
        XP2X = X(VarP2X);
        XP2Y = X(VarP2Y);
        
        XU1X = EST.XG(2*Nod1-1);
        XU1Y = EST.XG(2*Nod1);
        XU2X = EST.XG(2*Nod2-1);
        XU2Y = EST.XG(2*Nod2);
        
        X1 = EST.MatNod1(Nod1,1);
        Y1 = EST.MatNod1(Nod1,2);
        X2 = EST.MatNod1(Nod2,1);
        Y2 = EST.MatNod1(Nod2,2);
        
        DX = X2 - X1;
        DY = Y2 - Y1;
        
        DXX = DX*DX;
        DXY = DX*DY;
        DYY = DY*DY;
        
        LQ = DXX + DYY;
        LE = sqrt(LQ);
        LI = CD1/LE;
        LogLE = log(LE);
        
        R1X = X1 - XF;
        R1Y = Y1 - YF;
        R2X = X2 - XF;
        R2Y = Y2 - YF;
        
        Ge11 = CD0;  Ge12 = CD0;
        Ge13 = CD0;  Ge14 = CD0;
        Ge21 = CD0;  Ge22 = CD0;
        Ge23 = CD0;  Ge24 = CD0;
        
        He11 = CD0;  He12 = CD0;
        He13 = CD0;  He14 = CD0;
        He21 = CD0;  He22 = CD0;
        He23 = CD0;  He24 = CD0;
        
        Ra1 = sqrt(R1X*R1X + R1Y*R1Y);
        Ra2 = sqrt(R2X*R2X + R2Y*R2Y);
        
        if (Ra1 < TolRa)
            
            C4L = CM4*LI;
            
            C4DXX = C4L*DXX;
            C4DXY = C4L*DXY;
            C4DYY = C4L*DYY;
            
            CO1 = LogLE - CD6;
            CO2 = LogLE - CD5;
            
            C6LC1 = CM6*LE*CO1;
            C6LC2 = CM6*LE*CO2;
            
            Ge11 = C6LC1 - C4DXX;
            Ge22 = C6LC1 - C4DYY;
            Ge13 = C6LC2 - C4DXX;
            Ge24 = C6LC2 - C4DYY;
            Ge12 = -C4DXY;
            Ge14 = Ge12;
            Ge21 = Ge12;
            Ge23 = Ge12;
            
            He14 = CM7;
            He23 = -CM7;
            
        elseif (Ra2 < TolRa)
            
            C4L = CM4*LI;
            
            C4DXX = C4L*DXX;
            C4DXY = C4L*DXY;
            C4DYY = C4L*DYY;
            
            CO1 = LogLE - CD5;
            CO2 = LogLE - CD6;
            
            C6LC1 = CM6*LE*CO1;
            C6LC2 = CM6*LE*CO2;
            
            Ge11 = C6LC1 - C4DXX;
            Ge22 = C6LC1 - C4DYY;
            Ge13 = C6LC2 - C4DXX;
            Ge24 = C6LC2 - C4DYY;
            Ge12 = -C4DXY;
            Ge14 = Ge12;
            Ge21 = Ge12;
            Ge23 = Ge12;
            
            He12 = -CM7;
            He21 = CM7;
            
        else
            
            for PI = 1:NG
                
                N1 = CD5*(CD1 - XI(PI));
                N2 = CD5*(CD1 + XI(PI));
                
                RX = N1*R1X + N2*R2X;
                RY = N1*R1Y + N2*R2Y;
                
                RQ = RX*RX + RY*RY;
                RE = sqrt(RQ);
                RI = CD1/RE;
                LogRE = log(RE);
                
                dRX = RI*RX;
                dRY = RI*RY;
                
                dRN = dRX*DY - dRY*DX;
                dRT = dRX*DX + dRY*DY;
                
                C6R = CM5*LogRE;
                
                dRXX = dRX*dRX;
                dRXY = dRX*dRY;
                dRYY = dRY*dRY;
                
                UXX = C6R - CM3*dRXX;
                UYY = C6R - CM3*dRYY;
                UXY = -CM3*dRXY;
                
                C4XX = CM7 + CM2*dRXX;
                C4YY = CM7 + CM2*dRYY;
                C4XY = CM2*dRXY;
                
                dRNRI = RI*dRN;
                dRTRI = RI*dRT;
                
                PXX = C4XX*dRNRI;
                PYY = C4YY*dRNRI;
                
                PXY = C4XY*dRNRI + CM7*dRTRI;
                PYX = C4XY*dRNRI - CM7*dRTRI;
                
                W5 = CD5*WI(PI);
                
                WUXX = W5*LE*UXX;
                WUXY = W5*LE*UXY;
                WUYY = W5*LE*UYY;
                
                WPXX = W5*PXX;
                WPXY = W5*PXY;
                WPYX = W5*PYX;
                WPYY = W5*PYY;
                
                He11 = He11 + WPXX*N1;
                He12 = He12 + WPXY*N1;
                He13 = He13 + WPXX*N2;
                He14 = He14 + WPXY*N2;
                
                He21 = He21 + WPYX*N1;
                He22 = He22 + WPYY*N1;
                He23 = He23 + WPYX*N2;
                He24 = He24 + WPYY*N2;
                
                Ge11 = Ge11 + WUXX*N1;
                Ge12 = Ge12 + WUXY*N1;
                Ge13 = Ge13 + WUXX*N2;
                Ge14 = Ge14 + WUXY*N2;
                
                Ge21 = Ge21 + WUXY*N1;
                Ge22 = Ge22 + WUYY*N1;
                Ge23 = Ge23 + WUXY*N2;
                Ge24 = Ge24 + WUYY*N2;
                
            end
        end
        
        % Monto matrices MatH y MatG --------------------------------------
        
        % MatH(GrF,GrU) = MatH(GrF,GrU) + He;
        % MatG(GrF,GrP) = MatG(GrF,GrP) + Ge;
        
        G(ResF1) = G(ResF1) - Ge11*XP1X - Ge12*XP1Y - Ge13*XP2X - Ge14*XP2Y;
        G(ResF2) = G(ResF2) - Ge21*XP1X - Ge22*XP1Y - Ge23*XP2X - Ge24*XP2Y;
        G(ResF1) = G(ResF1) + He11*XU1X + He12*XU1Y + He13*XU2X + He14*XU2Y;
        G(ResF2) = G(ResF2) + He21*XU1X + He22*XU1Y + He23*XU2X + He24*XU2Y;
        
        % Desplazamientos rígidos -----------------------------------------
        
        % MatH(GrF,GrF) = MatH(GrF,GrF) - He(1:2,1:2) - He(1:2,3:4);
        
        G(ResF1) = G(ResF1) - He11*XUFX - He12*XUFY - He13*XUFX - He14*XUFY;
        G(ResF2) = G(ResF2) - He21*XUFX - He22*XUFY - He23*XUFX - He24*XUFY;
        
    end
    
    for Ele = 1:EST.NumEle2
        
        Nod1 = EST.MatEle2(Ele,1);
        Nod2 = EST.MatEle2(Ele,2);
        
        VarU1X = MatXU(2*Nod1-1);
        VarU1Y = MatXU(2*Nod1);
        VarU2X = MatXU(2*Nod2-1);
        VarU2Y = MatXU(2*Nod2);
        
        XU1X = X(VarU1X);
        XU1Y = X(VarU1Y);
        XU2X = X(VarU2X);
        XU2Y = X(VarU2Y);
        
        VarT1 = MatXT(Nod1);
        VarT2 = MatXT(Nod2);
        
        VN1X = EST.VZ(Nod1,1);
        VN1Y = EST.VZ(Nod1,2);
        VN2X = EST.VZ(Nod2,1);
        VN2Y = EST.VZ(Nod2,2);
        
        X1 = EST.MatNod2(Nod1,1) + VN1X*X(VarT1);
        Y1 = EST.MatNod2(Nod1,2) + VN1Y*X(VarT1);
        X2 = EST.MatNod2(Nod2,1) + VN2X*X(VarT2);
        Y2 = EST.MatNod2(Nod2,2) + VN2Y*X(VarT2);
        
        DX = X2 - X1;
        DY = Y2 - Y1;
        
        DXX = DX*DX;
        DXY = DX*DY;
        DYY = DY*DY;
        
        LQ = DXX + DYY;
        LE = sqrt(LQ);
        LI = CD1/LE;
        LogLE = log(LE);
        
        % Presion exterior ------------------------------------------------
        XP1 = -EST.p0; % presion impuesta
        XP2 = -EST.p0; % presion impuesta
        
        XP1 = XP1 + EST.bA*Y1; % resto solucion particular
        XP2 = XP2 + EST.bA*Y2; % resto solucion particular
        
        % Terminos de la curvatura
        EleA = Ele-1; if (EleA < 1), EleA = EST.NumEle2; end
        EleB = Ele+1; if (EleB > EST.NumEle2), EleB = 1; end
        
        NodA = EST.MatEle2(EleA,1);
        NodB = EST.MatEle2(EleB,2);
        
        VarTA = MatXT(NodA);
        VarTB = MatXT(NodB);
        
        VNAX = EST.VZ(NodA,1);
        VNAY = EST.VZ(NodA,2);
        VNBX = EST.VZ(NodB,1);
        VNBY = EST.VZ(NodB,2);
        
        XA = EST.MatNod2(NodA,1) + VNAX*X(VarTA);
        YA = EST.MatNod2(NodA,2) + VNAY*X(VarTA);
        XB = EST.MatNod2(NodB,1) + VNBX*X(VarTB);
        YB = EST.MatNod2(NodB,2) + VNBY*X(VarTB);
        
        % Longitudes
        DXA = X1 - XA; DXA3 = X2 - XA;
        DYA = Y1 - YA; DYA3 = Y2 - YA;
        DXB = XB - X2; DXB3 = XB - X1;
        DYB = YB - Y2; DYB3 = YB - Y1;
        
        LQA = DXA*DXA + DYA*DYA;
        LQA3 = DXA3*DXA3 + DYA3*DYA3;
        LQB = DXB*DXB + DYB*DYB;
        LQB3 = DXB3*DXB3 + DYB3*DYB3;
        
        LEA = sqrt(LQA);
        LIA = CD1/LEA;
        LEA3 = sqrt(LQA3);
        LIA3 = CD1/LEA3;
        LEB = sqrt(LQB);
        LIB = CD1/LEB;
        LEB3 = sqrt(LQB3);
        LIB3 = CD1/LEB3;
        
        SGA = sign(DXA*DY - DYA*DX);
        SGB = sign(DX*DYB - DY*DXB);
        
        LA1 = LEA + LE + LEA3;   LA2 = LEA + LE - LEA3;
        LA3 = LEA3 + LEA - LE;   LA4 = LE + LEA3 - LEA;
        
        LA12 = LA1*LA2;          LA34 = LA3*LA4;
        LA14 = LA12*LA34;        LIAS = LI*LIA;
        LAR4 = sqrt(LA14);       LISA = LIAS*LIA3;
        
        C1 = LAR4*LISA; % Curvatura nodo 1
        C1 = SGA*C1;    % Curvatura nodo 1
        
        LB1 = LEB + LE + LEB3;   LB2 = LEB + LE - LEB3;
        LB3 = LEB3 + LEB - LE;   LB4 = LE + LEB3 - LEB;
        
        LB12 = LB1*LB2;          LB34 = LB3*LB4;
        LB14 = LB12*LB34;        LIBS = LI*LIB;
        LBR4 = sqrt(LB14);       LISB = LIBS*LIB3;
        
        C2 = LBR4*LISB; % Curvatura nodo 2
        C2 = SGB*C2;    % Curvatura nodo 2
        
        XP1 = XP1 - EST.Gamma*C1;
        XP2 = XP2 - EST.Gamma*C2;
        
        NX = +LI*DY;
        NY = -LI*DX;
        
        XP1X = XP1*NX;
        XP1Y = XP1*NY;
        XP2X = XP2*NX;
        XP2Y = XP2*NY;
        % -----------------------------------------------------------------
        
        R1X = X1 - XF;
        R1Y = Y1 - YF;
        R2X = X2 - XF;
        R2Y = Y2 - YF;
        
        Ge11 = CD0;  Ge12 = CD0;
        Ge13 = CD0;  Ge14 = CD0;
        Ge21 = CD0;  Ge22 = CD0;
        Ge23 = CD0;  Ge24 = CD0;
        
        He11 = CD0;  He12 = CD0;
        He13 = CD0;  He14 = CD0;
        He21 = CD0;  He22 = CD0;
        He23 = CD0;  He24 = CD0;
        
        Ra1 = sqrt(R1X*R1X + R1Y*R1Y);
        Ra2 = sqrt(R2X*R2X + R2Y*R2Y);
        
        if (Ra1 < TolRa)
            
            C4L = CM4*LI;
            
            C4DXX = C4L*DXX;
            C4DXY = C4L*DXY;
            C4DYY = C4L*DYY;
            
            CO1 = LogLE - CD6;
            CO2 = LogLE - CD5;
            
            C6LC1 = CM6*LE*CO1;
            C6LC2 = CM6*LE*CO2;
            
            Ge11 = C6LC1 - C4DXX;
            Ge22 = C6LC1 - C4DYY;
            Ge13 = C6LC2 - C4DXX;
            Ge24 = C6LC2 - C4DYY;
            Ge12 = -C4DXY;
            Ge14 = Ge12;
            Ge21 = Ge12;
            Ge23 = Ge12;
            
            He14 = CM7;
            He23 = -CM7;
            
        elseif (Ra2 < TolRa)
            
            C4L = CM4*LI;
            
            C4DXX = C4L*DXX;
            C4DXY = C4L*DXY;
            C4DYY = C4L*DYY;
            
            CO1 = LogLE - CD5;
            CO2 = LogLE - CD6;
            
            C6LC1 = CM6*LE*CO1;
            C6LC2 = CM6*LE*CO2;
            
            Ge11 = C6LC1 - C4DXX;
            Ge22 = C6LC1 - C4DYY;
            Ge13 = C6LC2 - C4DXX;
            Ge24 = C6LC2 - C4DYY;
            Ge12 = -C4DXY;
            Ge14 = Ge12;
            Ge21 = Ge12;
            Ge23 = Ge12;
            
            He12 = -CM7;
            He21 = CM7;
            
        else
            
            for PI = 1:NG
                
                N1 = CD5*(CD1 - XI(PI));
                N2 = CD5*(CD1 + XI(PI));
                
                RX = N1*R1X + N2*R2X;
                RY = N1*R1Y + N2*R2Y;
                
                RQ = RX*RX + RY*RY;
                RE = sqrt(RQ);
                RI = CD1/RE;
                LogRE = log(RE);
                
                dRX = RI*RX;
                dRY = RI*RY;
                
                dRN = dRX*DY - dRY*DX;
                dRT = dRX*DX + dRY*DY;
                
                C6R = CM5*LogRE;
                
                dRXX = dRX*dRX;
                dRXY = dRX*dRY;
                dRYY = dRY*dRY;
                
                UXX = C6R - CM3*dRXX;
                UYY = C6R - CM3*dRYY;
                UXY = -CM3*dRXY;
                
                C4XX = CM7 + CM2*dRXX;
                C4YY = CM7 + CM2*dRYY;
                C4XY = CM2*dRXY;
                
                dRNRI = RI*dRN;
                dRTRI = RI*dRT;
                
                PXX = C4XX*dRNRI;
                PYY = C4YY*dRNRI;
                
                PXY = C4XY*dRNRI + CM7*dRTRI;
                PYX = C4XY*dRNRI - CM7*dRTRI;
                
                W5 = CD5*WI(PI);
                
                WUXX = W5*LE*UXX;
                WUXY = W5*LE*UXY;
                WUYY = W5*LE*UYY;
                
                WPXX = W5*PXX;
                WPXY = W5*PXY;
                WPYX = W5*PYX;
                WPYY = W5*PYY;
                
                He11 = He11 + WPXX*N1;
                He12 = He12 + WPXY*N1;
                He13 = He13 + WPXX*N2;
                He14 = He14 + WPXY*N2;
                
                He21 = He21 + WPYX*N1;
                He22 = He22 + WPYY*N1;
                He23 = He23 + WPYX*N2;
                He24 = He24 + WPYY*N2;
                
                Ge11 = Ge11 + WUXX*N1;
                Ge12 = Ge12 + WUXY*N1;
                Ge13 = Ge13 + WUXX*N2;
                Ge14 = Ge14 + WUXY*N2;
                
                Ge21 = Ge21 + WUXY*N1;
                Ge22 = Ge22 + WUYY*N1;
                Ge23 = Ge23 + WUXY*N2;
                Ge24 = Ge24 + WUYY*N2;
                
            end
        end
        
        % Monto matrices MatH y MatG --------------------------------------
        
        % MatH(GrF,GrU) = MatH(GrF,GrU) + He;
        % MatG(GrF,GrP) = MatG(GrF,GrP) + Ge;
        
        G(ResF1) = G(ResF1) - Ge11*XP1X - Ge12*XP1Y - Ge13*XP2X - Ge14*XP2Y;
        G(ResF2) = G(ResF2) - Ge21*XP1X - Ge22*XP1Y - Ge23*XP2X - Ge24*XP2Y;
        G(ResF1) = G(ResF1) + He11*XU1X + He12*XU1Y + He13*XU2X + He14*XU2Y;
        G(ResF2) = G(ResF2) + He21*XU1X + He22*XU1Y + He23*XU2X + He24*XU2Y;
        
        % Desplazamientos rígidos -----------------------------------------
        
        % MatH(GrF,GrF) = MatH(GrF,GrF) - He(1:2,1:2) - He(1:2,3:4);
        
        G(ResF1) = G(ResF1) - He11*XUFX - He12*XUFY - He13*XUFX - He14*XUFY;
        G(ResF2) = G(ResF2) - He21*XUFX - He22*XUFY - He23*XUFX - He24*XUFY;
        
    end
    
    % Variable lamda en el sistema ----------------------------------------
    VarLAM = MatXL;
    
    XLAM = X(VarLAM);
    
    G(ResF1) = G(ResF1) + EST.PN(2*NodF-1)*XLAM;
    G(ResF2) = G(ResF2) + EST.PN(2*NodF)*XLAM;
    
    VarPFX = MatXP(2*NodF-1);
    VarPFY = MatXP(2*NodF);
    
    XPFX = X(VarPFX);
    XPFY = X(VarPFY);
    
    G(ResLAM) = G(ResLAM) + EST.PN(2*NodF-1)*XPFX;
    G(ResLAM) = G(ResLAM) + EST.PN(2*NodF)*XPFY;
    
end

for NodF = 1:EST.NumNod2
    
    VarTF = MatXT(NodF);
    
    VNFX = EST.VZ(NodF,1);
    VNFY = EST.VZ(NodF,2);
    
    XF = EST.MatNod2(NodF,1) + VNFX*X(VarTF);
    YF = EST.MatNod2(NodF,2) + VNFY*X(VarTF);
    
    VarUFX = MatXU(2*NodF-1);
    VarUFY = MatXU(2*NodF);
    
    XUFX = X(VarUFX);
    XUFY = X(VarUFY);
    
    ResF1 = 2*EST.NumNod1 + 2*NodF-1;
    ResF2 = 2*EST.NumNod1 + 2*NodF;
    
    for Ele = 1:EST.NumEle1
        
        Nod1 = EST.MatEle1(Ele,1);
        Nod2 = EST.MatEle1(Ele,2);
        
        VarP1X = MatXP(2*Nod1-1);
        VarP1Y = MatXP(2*Nod1);
        VarP2X = MatXP(2*Nod2-1);
        VarP2Y = MatXP(2*Nod2);
        
        XP1X = X(VarP1X);
        XP1Y = X(VarP1Y);
        XP2X = X(VarP2X);
        XP2Y = X(VarP2Y);
        
        XU1X = EST.XG(2*Nod1-1);
        XU1Y = EST.XG(2*Nod1);
        XU2X = EST.XG(2*Nod2-1);
        XU2Y = EST.XG(2*Nod2);
        
        X1 = EST.MatNod1(Nod1,1);
        Y1 = EST.MatNod1(Nod1,2);
        X2 = EST.MatNod1(Nod2,1);
        Y2 = EST.MatNod1(Nod2,2);
        
        DX = X2 - X1;
        DY = Y2 - Y1;
        
        DXX = DX*DX;
        DXY = DX*DY;
        DYY = DY*DY;
        
        LQ = DXX + DYY;
        LE = sqrt(LQ);
        LI = CD1/LE;
        LogLE = log(LE);
        
        R1X = X1 - XF;
        R1Y = Y1 - YF;
        R2X = X2 - XF;
        R2Y = Y2 - YF;
        
        Ge11 = CD0;  Ge12 = CD0;
        Ge13 = CD0;  Ge14 = CD0;
        Ge21 = CD0;  Ge22 = CD0;
        Ge23 = CD0;  Ge24 = CD0;
        
        He11 = CD0;  He12 = CD0;
        He13 = CD0;  He14 = CD0;
        He21 = CD0;  He22 = CD0;
        He23 = CD0;  He24 = CD0;
        
        Ra1 = sqrt(R1X*R1X + R1Y*R1Y);
        Ra2 = sqrt(R2X*R2X + R2Y*R2Y);
        
        if (Ra1 < TolRa)
            
            C4L = CM4*LI;
            
            C4DXX = C4L*DXX;
            C4DXY = C4L*DXY;
            C4DYY = C4L*DYY;
            CO1 = LogLE - CD6;
            CO2 = LogLE - CD5;
            
            C6LC1 = CM6*LE*CO1;
            C6LC2 = CM6*LE*CO2;
            
            Ge11 = C6LC1 - C4DXX;
            Ge22 = C6LC1 - C4DYY;
            Ge13 = C6LC2 - C4DXX;
            Ge24 = C6LC2 - C4DYY;
            Ge12 = -C4DXY;
            Ge14 = Ge12;
            Ge21 = Ge12;
            Ge23 = Ge12;
            
            He14 = CM7;
            He23 = -CM7;
            
        elseif (Ra2 < TolRa)
            
            C4L = CM4*LI;
            
            C4DXX = C4L*DXX;
            C4DXY = C4L*DXY;
            C4DYY = C4L*DYY;
            
            CO1 = LogLE - CD5;
            CO2 = LogLE - CD6;
            
            C6LC1 = CM6*LE*CO1;
            C6LC2 = CM6*LE*CO2;
            
            Ge11 = C6LC1 - C4DXX;
            Ge22 = C6LC1 - C4DYY;
            Ge13 = C6LC2 - C4DXX;
            Ge24 = C6LC2 - C4DYY;
            Ge12 = -C4DXY;
            Ge14 = Ge12;
            Ge21 = Ge12;
            Ge23 = Ge12;
            
            He12 = -CM7;
            He21 = CM7;
            
        else
            
            for PI = 1:NG
                
                N1 = CD5*(CD1 - XI(PI));
                N2 = CD5*(CD1 + XI(PI));
                
                RX = N1*R1X + N2*R2X;
                RY = N1*R1Y + N2*R2Y;
                
                RQ = RX*RX + RY*RY;
                RE = sqrt(RQ);
                RI = CD1/RE;
                LogRE = log(RE);
                
                dRX = RI*RX;
                dRY = RI*RY;
                
                dRN = dRX*DY - dRY*DX;
                dRT = dRX*DX + dRY*DY;
                
                C6R = CM5*LogRE;
                
                dRXX = dRX*dRX;
                dRXY = dRX*dRY;
                dRYY = dRY*dRY;
                
                UXX = C6R - CM3*dRXX;
                UYY = C6R - CM3*dRYY;
                UXY = -CM3*dRXY;
                
                C4XX = CM7 + CM2*dRXX;
                C4YY = CM7 + CM2*dRYY;
                C4XY = CM2*dRXY;
                
                dRNRI = RI*dRN;
                dRTRI = RI*dRT;
                
                PXX = C4XX*dRNRI;
                PYY = C4YY*dRNRI;
                
                PXY = C4XY*dRNRI + CM7*dRTRI;
                PYX = C4XY*dRNRI - CM7*dRTRI;
                
                W5 = CD5*WI(PI);
                
                WUXX = W5*LE*UXX;
                WUXY = W5*LE*UXY;
                WUYY = W5*LE*UYY;
                
                WPXX = W5*PXX;
                WPXY = W5*PXY;
                WPYX = W5*PYX;
                WPYY = W5*PYY;
                
                He11 = He11 + WPXX*N1;
                He12 = He12 + WPXY*N1;
                He13 = He13 + WPXX*N2;
                He14 = He14 + WPXY*N2;
                
                He21 = He21 + WPYX*N1;
                He22 = He22 + WPYY*N1;
                He23 = He23 + WPYX*N2;
                He24 = He24 + WPYY*N2;
                
                Ge11 = Ge11 + WUXX*N1;
                Ge12 = Ge12 + WUXY*N1;
                Ge13 = Ge13 + WUXX*N2;
                Ge14 = Ge14 + WUXY*N2;
                
                Ge21 = Ge21 + WUXY*N1;
                Ge22 = Ge22 + WUYY*N1;
                Ge23 = Ge23 + WUXY*N2;
                Ge24 = Ge24 + WUYY*N2;
                
            end
        end
        
        % Monto matrices MatH y MatG --------------------------------------
        
        % MatH(GrF,GrU) = MatH(GrF,GrU) + He;
        % MatG(GrF,GrP) = MatG(GrF,GrP) + Ge;
        
        G(ResF1) = G(ResF1) - Ge11*XP1X - Ge12*XP1Y - Ge13*XP2X - Ge14*XP2Y;
        G(ResF2) = G(ResF2) - Ge21*XP1X - Ge22*XP1Y - Ge23*XP2X - Ge24*XP2Y;
        G(ResF1) = G(ResF1) + He11*XU1X + He12*XU1Y + He13*XU2X + He14*XU2Y;
        G(ResF2) = G(ResF2) + He21*XU1X + He22*XU1Y + He23*XU2X + He24*XU2Y;
        
        % Desplazamientos rígidos -----------------------------------------
        
        % MatH(GrF,GrF) = MatH(GrF,GrF) - He(1:2,1:2) - He(1:2,3:4);
        
        G(ResF1) = G(ResF1) - He11*XUFX - He12*XUFY - He13*XUFX - He14*XUFY;
        G(ResF2) = G(ResF2) - He21*XUFX - He22*XUFY - He23*XUFX - He24*XUFY;
        
    end
    
    for Ele = 1:EST.NumEle2
        
        Nod1 = EST.MatEle2(Ele,1);
        Nod2 = EST.MatEle2(Ele,2);
        
        VarU1X = MatXU(2*Nod1-1);
        VarU1Y = MatXU(2*Nod1);
        VarU2X = MatXU(2*Nod2-1);
        VarU2Y = MatXU(2*Nod2);
        
        XU1X = X(VarU1X);
        XU1Y = X(VarU1Y);
        XU2X = X(VarU2X);
        XU2Y = X(VarU2Y);
        
        VarT1 = MatXT(Nod1);
        VarT2 = MatXT(Nod2);
        
        VN1X = EST.VZ(Nod1,1);
        VN1Y = EST.VZ(Nod1,2);
        VN2X = EST.VZ(Nod2,1);
        VN2Y = EST.VZ(Nod2,2);
        
        X1 = EST.MatNod2(Nod1,1) + VN1X*X(VarT1);
        Y1 = EST.MatNod2(Nod1,2) + VN1Y*X(VarT1);
        X2 = EST.MatNod2(Nod2,1) + VN2X*X(VarT2);
        Y2 = EST.MatNod2(Nod2,2) + VN2Y*X(VarT2);
        
        DX = X2 - X1;
        DY = Y2 - Y1;
        
        DXX = DX*DX;
        DXY = DX*DY;
        DYY = DY*DY;
        
        LQ = DXX + DYY;
        LE = sqrt(LQ);
        LI = CD1/LE;
        LogLE = log(LE);
        
        % Presion exterior ------------------------------------------------
        XP1 = -EST.p0; % presion impuesta
        XP2 = -EST.p0; % presion impuesta
        
        XP1 = XP1 + EST.bA*Y1; % resto solucion particular
        XP2 = XP2 + EST.bA*Y2; % resto solucion particular
        
        % Terminos de la curvatura
        EleA = Ele-1; if (EleA < 1), EleA = EST.NumEle2; end
        EleB = Ele+1; if (EleB > EST.NumEle2), EleB = 1; end
        
        NodA = EST.MatEle2(EleA,1);
        NodB = EST.MatEle2(EleB,2);
        
        VarTA = MatXT(NodA);
        VarTB = MatXT(NodB);
        
        VNAX = EST.VZ(NodA,1);
        VNAY = EST.VZ(NodA,2);
        VNBX = EST.VZ(NodB,1);
        VNBY = EST.VZ(NodB,2);
        
        XA = EST.MatNod2(NodA,1) + VNAX*X(VarTA);
        YA = EST.MatNod2(NodA,2) + VNAY*X(VarTA);
        XB = EST.MatNod2(NodB,1) + VNBX*X(VarTB);
        YB = EST.MatNod2(NodB,2) + VNBY*X(VarTB);
        
        % Longitudes
        DXA = X1 - XA; DXA3 = X2 - XA;
        DYA = Y1 - YA; DYA3 = Y2 - YA;
        DXB = XB - X2; DXB3 = XB - X1;
        DYB = YB - Y2; DYB3 = YB - Y1;
        
        LQA = DXA*DXA + DYA*DYA;
        LQA3 = DXA3*DXA3 + DYA3*DYA3;
        LQB = DXB*DXB + DYB*DYB;
        LQB3 = DXB3*DXB3 + DYB3*DYB3;
        
        LEA = sqrt(LQA);
        LIA = CD1/LEA;
        LEA3 = sqrt(LQA3);
        LIA3 = CD1/LEA3;
        LEB = sqrt(LQB);
        LIB = CD1/LEB;
        LEB3 = sqrt(LQB3);
        LIB3 = CD1/LEB3;
        
        SGA = sign(DXA*DY - DYA*DX);
        SGB = sign(DX*DYB - DY*DXB);
        
        LA1 = LEA + LE + LEA3;   LA2 = LEA + LE - LEA3;
        LA3 = LEA3 + LEA - LE;   LA4 = LE + LEA3 - LEA;
        
        LA12 = LA1*LA2;          LA34 = LA3*LA4;
        LA14 = LA12*LA34;        LIAS = LI*LIA;
        LAR4 = sqrt(LA14);       LISA = LIAS*LIA3;
        
        C1 = LAR4*LISA; % Curvatura nodo 1
        C1 = SGA*C1;    % Curvatura nodo 1
        
        LB1 = LEB + LE + LEB3;   LB2 = LEB + LE - LEB3;
        LB3 = LEB3 + LEB - LE;   LB4 = LE + LEB3 - LEB;
        
        LB12 = LB1*LB2;          LB34 = LB3*LB4;
        LB14 = LB12*LB34;        LIBS = LI*LIB;
        LBR4 = sqrt(LB14);       LISB = LIBS*LIB3;
        
        C2 = LBR4*LISB; % Curvatura nodo 2
        C2 = SGB*C2;    % Curvatura nodo 2
        
        XP1 = XP1 - EST.Gamma*C1;
        XP2 = XP2 - EST.Gamma*C2;
        
        NX = +LI*DY;
        NY = -LI*DX;
        
        XP1X = XP1*NX;
        XP1Y = XP1*NY;
        XP2X = XP2*NX;
        XP2Y = XP2*NY;
        % -----------------------------------------------------------------
        
        R1X = X1 - XF;
        R1Y = Y1 - YF;
        R2X = X2 - XF;
        R2Y = Y2 - YF;
        
        Ge11 = CD0;  Ge12 = CD0;
        Ge13 = CD0;  Ge14 = CD0;
        Ge21 = CD0;  Ge22 = CD0;
        Ge23 = CD0;  Ge24 = CD0;
        
        He11 = CD0;  He12 = CD0;
        He13 = CD0;  He14 = CD0;
        He21 = CD0;  He22 = CD0;
        He23 = CD0;  He24 = CD0;
        
        Ra1 = sqrt(R1X*R1X + R1Y*R1Y);
        Ra2 = sqrt(R2X*R2X + R2Y*R2Y);
        
        if (Ra1 < TolRa)
            
            C4L = CM4*LI;
            
            C4DXX = C4L*DXX;
            C4DXY = C4L*DXY;
            C4DYY = C4L*DYY;
            CO1 = LogLE - CD6;
            CO2 = LogLE - CD5;
            
            C6LC1 = CM6*LE*CO1;
            C6LC2 = CM6*LE*CO2;
            
            Ge11 = C6LC1 - C4DXX;
            Ge22 = C6LC1 - C4DYY;
            Ge13 = C6LC2 - C4DXX;
            Ge24 = C6LC2 - C4DYY;
            Ge12 = -C4DXY;
            Ge14 = Ge12;
            Ge21 = Ge12;
            Ge23 = Ge12;
            
            He14 = CM7;
            He23 = -CM7;
            
        elseif (Ra2 < TolRa)
            
            C4L = CM4*LI;
            
            C4DXX = C4L*DXX;
            C4DXY = C4L*DXY;
            C4DYY = C4L*DYY;
            
            CO1 = LogLE - CD5;
            CO2 = LogLE - CD6;
            
            C6LC1 = CM6*LE*CO1;
            C6LC2 = CM6*LE*CO2;
            
            Ge11 = C6LC1 - C4DXX;
            Ge22 = C6LC1 - C4DYY;
            Ge13 = C6LC2 - C4DXX;
            Ge24 = C6LC2 - C4DYY;
            Ge12 = -C4DXY;
            Ge14 = Ge12;
            Ge21 = Ge12;
            Ge23 = Ge12;
            
            He12 = -CM7;
            He21 = CM7;
            
        else
            
            for PI = 1:NG
                
                N1 = CD5*(CD1 - XI(PI));
                N2 = CD5*(CD1 + XI(PI));
                
                RX = N1*R1X + N2*R2X;
                RY = N1*R1Y + N2*R2Y;
                
                RQ = RX*RX + RY*RY;
                RE = sqrt(RQ);
                RI = CD1/RE;
                LogRE = log(RE);
                
                dRX = RI*RX;
                dRY = RI*RY;
                
                dRN = dRX*DY - dRY*DX;
                dRT = dRX*DX + dRY*DY;
                
                C6R = CM5*LogRE;
                
                dRXX = dRX*dRX;
                dRXY = dRX*dRY;
                dRYY = dRY*dRY;
                
                UXX = C6R - CM3*dRXX;
                UYY = C6R - CM3*dRYY;
                UXY = -CM3*dRXY;
                
                C4XX = CM7 + CM2*dRXX;
                C4YY = CM7 + CM2*dRYY;
                C4XY = CM2*dRXY;
                
                dRNRI = RI*dRN;
                dRTRI = RI*dRT;
                
                PXX = C4XX*dRNRI;
                PYY = C4YY*dRNRI;
                
                PXY = C4XY*dRNRI + CM7*dRTRI;
                PYX = C4XY*dRNRI - CM7*dRTRI;
                
                W5 = CD5*WI(PI);
                
                WUXX = W5*LE*UXX;
                WUXY = W5*LE*UXY;
                WUYY = W5*LE*UYY;
                
                WPXX = W5*PXX;
                WPXY = W5*PXY;
                WPYX = W5*PYX;
                WPYY = W5*PYY;
                
                He11 = He11 + WPXX*N1;
                He12 = He12 + WPXY*N1;
                He13 = He13 + WPXX*N2;
                He14 = He14 + WPXY*N2;
                
                He21 = He21 + WPYX*N1;
                He22 = He22 + WPYY*N1;
                He23 = He23 + WPYX*N2;
                He24 = He24 + WPYY*N2;
                
                Ge11 = Ge11 + WUXX*N1;
                Ge12 = Ge12 + WUXY*N1;
                Ge13 = Ge13 + WUXX*N2;
                Ge14 = Ge14 + WUXY*N2;
                
                Ge21 = Ge21 + WUXY*N1;
                Ge22 = Ge22 + WUYY*N1;
                Ge23 = Ge23 + WUXY*N2;
                Ge24 = Ge24 + WUYY*N2;
                
            end
        end
        
        % Monto matrices MatH y MatG --------------------------------------
        
        % MatH(GrF,GrU) = MatH(GrF,GrU) + He;
        % MatG(GrF,GrP) = MatG(GrF,GrP) + Ge;
        
        G(ResF1) = G(ResF1) - Ge11*XP1X - Ge12*XP1Y - Ge13*XP2X - Ge14*XP2Y;
        G(ResF2) = G(ResF2) - Ge21*XP1X - Ge22*XP1Y - Ge23*XP2X - Ge24*XP2Y;
        G(ResF1) = G(ResF1) + He11*XU1X + He12*XU1Y + He13*XU2X + He14*XU2Y;
        G(ResF2) = G(ResF2) + He21*XU1X + He22*XU1Y + He23*XU2X + He24*XU2Y;
        
        % Desplazamientos rígidos -----------------------------------------
        
        % MatH(GrF,GrF) = MatH(GrF,GrF) - He(1:2,1:2) - He(1:2,3:4);
        
        G(ResF1) = G(ResF1) - He11*XUFX - He12*XUFY - He13*XUFX - He14*XUFY;
        G(ResF2) = G(ResF2) - He21*XUFX - He22*XUFY - He23*XUFX - He24*XUFY;
        
    end
end

for Ele = 1:EST.NumNod2
    
    Nod1 = EST.MatEle2(Ele,1);
    Nod2 = EST.MatEle2(Ele,2);
    
    ResF1 = 2*EST.NumNod1 + 2*EST.NumNod2 + Nod1;
    ResF2 = 2*EST.NumNod1 + 2*EST.NumNod2 + Nod2;
    
    VarT1 = MatXT(Nod1);
    VarT2 = MatXT(Nod2);
    
    VN1X = EST.VZ(Nod1,1);
    VN1Y = EST.VZ(Nod1,2);
    VN2X = EST.VZ(Nod2,1);
    VN2Y = EST.VZ(Nod2,2);
    
    X1 = EST.MatNod2(Nod1,1) + VN1X*X(VarT1);
    Y1 = EST.MatNod2(Nod1,2) + VN1Y*X(VarT1);
    X2 = EST.MatNod2(Nod2,1) + VN2X*X(VarT2);
    Y2 = EST.MatNod2(Nod2,2) + VN2Y*X(VarT2);
    
    VarU1X = MatXU(2*Nod1-1);
    VarU1Y = MatXU(2*Nod1);
    VarU2X = MatXU(2*Nod2-1);
    VarU2Y = MatXU(2*Nod2);
    
    % Solucion particular
    XU1XP = CD0;
    XU1YP = CD0;
    XU2XP = CD0;
    XU2YP = CD0;

    XU1X = X(VarU1X) + XU1XP; % sumo solucion particular
    XU1Y = X(VarU1Y) + XU1YP; % sumo solucion particular
    XU2X = X(VarU2X) + XU2XP; % sumo solucion particular
    XU2Y = X(VarU2Y) + XU2YP; % sumo solucion particular
    
    DX = X2 - X1;
    DY = Y2 - Y1;
    
    U1N = XU1X*DY - XU1Y*DX;
    U2N = XU2X*DY - XU2Y*DX;
    
    G(ResF1) = G(ResF1) + CI3*U1N + CI6*U2N; % integro UN*N1
    G(ResF2) = G(ResF2) + CI6*U1N + CI3*U2N; % integro UN*N2
    
end

end
