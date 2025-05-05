function gG = gfunFU(EST)

% Unmodified function derivative (grad F)

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
MatXL = 2*EST.NumNod1+3*EST.NumNod2+1;

ResLAM = EST.NumRes; % Restriccion de variable lambda

% G = zeros(EST.NumRes,1);

gG = zeros(EST.NumRes,EST.NumVar);

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
        
        % XP1X = X(VarP1X);
        % XP1Y = X(VarP1Y);
        % XP2X = X(VarP2X);
        % XP2Y = X(VarP2Y);
        
        gp1xXP1X = CD1;
        gp1yXP1Y = CD1;
        gp2xXP2X = CD1;
        gp2yXP2Y = CD1;
        
        % XU1X = EST.XG(2*Nod1-1);
        % XU1Y = EST.XG(2*Nod1);
        % XU2X = EST.XG(2*Nod2-1);
        % XU2Y = EST.XG(2*Nod2);
        
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
            
            % He14 = CM7;
            % He23 = -CM7;
            
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
            
            % He12 = -CM7;
            % He21 = CM7;
            
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
        
        % G(ResF1) = G(ResF1) - Ge11*XP1X - Ge12*XP1Y - Ge13*XP2X - Ge14*XP2Y;
        % G(ResF2) = G(ResF2) - Ge21*XP1X - Ge22*XP1Y - Ge23*XP2X - Ge24*XP2Y;
        % G(ResF1) = G(ResF1) + He11*XU1X + He12*XU1Y + He13*XU2X + He14*XU2Y;
        % G(ResF2) = G(ResF2) + He21*XU1X + He22*XU1Y + He23*XU2X + He24*XU2Y;
        
        gG(ResF1,VarP1X) = gG(ResF1,VarP1X) - Ge11*gp1xXP1X;
        gG(ResF1,VarP1Y) = gG(ResF1,VarP1Y) - Ge12*gp1yXP1Y;
        gG(ResF1,VarP2X) = gG(ResF1,VarP2X) - Ge13*gp2xXP2X;
        gG(ResF1,VarP2Y) = gG(ResF1,VarP2Y) - Ge14*gp2yXP2Y;
        gG(ResF2,VarP1X) = gG(ResF2,VarP1X) - Ge21*gp1xXP1X;
        gG(ResF2,VarP1Y) = gG(ResF2,VarP1Y) - Ge22*gp1yXP1Y;
        gG(ResF2,VarP2X) = gG(ResF2,VarP2X) - Ge23*gp2xXP2X;
        gG(ResF2,VarP2Y) = gG(ResF2,VarP2Y) - Ge24*gp2yXP2Y;
        
        % Desplazamientos rígidos -----------------------------------------
        
        % MatH(GrF,GrF) = MatH(GrF,GrF) - He(1:2,1:2) - He(1:2,3:4);
        
        % G(ResF1) = G(ResF1) - He11*XUFX - He12*XUFY - He13*XUFX - He14*XUFY;
        % G(ResF2) = G(ResF2) - He21*XUFX - He22*XUFY - He23*XUFX - He24*XUFY;
        
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
        
        gu1xXU1X = CD1;
        gu1yXU1Y = CD1;
        gu2xXU2X = CD1;
        gu2yXU2Y = CD1;
        
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
        
        gt1X1 = VN1X;
        gt1Y1 = VN1Y;
        gt2X2 = VN2X;
        gt2Y2 = VN2Y;
        
        DX = X2 - X1;
        DY = Y2 - Y1;
        
        gt1DX = -gt1X1;
        gt2DX =  gt2X2;
        gt1DY = -gt1Y1;
        gt2DY =  gt2Y2;
        
        DXX = DX*DX;
        DXY = DX*DY;
        DYY = DY*DY;
        
        gt1DXX = CD2*DX*gt1DX;
        gt2DXX = CD2*DX*gt2DX;
        gt1DXY = DX*gt1DY + DY*gt1DX;
        gt2DXY = DX*gt2DY + DY*gt2DX;
        gt1DYY = CD2*DY*gt1DY;
        gt2DYY = CD2*DY*gt2DY;
        
        LQ = DXX + DYY;
        LE = sqrt(LQ);
        LI = CD1/LE;
        LogLE = log(LE);
        
        gt1LQ = gt1DXX + gt1DYY;
        gt2LQ = gt2DXX + gt2DYY;
        gt1LE = CD5*LI*gt1LQ;
        gt2LE = CD5*LI*gt2LQ;
        gt1LI = -LI*LI*gt1LE;
        gt2LI = -LI*LI*gt2LE;
        gt1LogLE = LI*gt1LE;
        gt2LogLE = LI*gt2LE;
        
        % Presion exterior ------------------------------------------------
        XP1 = -EST.p0; % presion impuesta
        XP2 = -EST.p0; % presion impuesta
        
        gt1XP1 = CD0;
        gt2XP1 = CD0;
        gtaXP1 = CD0;
        gt1XP2 = CD0;
        gt2XP2 = CD0;
        gtbXP2 = CD0;
        
        XP1 = XP1 + EST.bA*Y1; % resto solucion particular
        XP2 = XP2 + EST.bA*Y2; % resto solucion particular
        
        gt1XP1 = gt1XP1 + EST.bA*gt1Y1;
        gt2XP2 = gt2XP2 + EST.bA*gt2Y2;
        
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
        
        gtaXA = VNAX;
        gtaYA = VNAY;
        gtbXB = VNBX;
        gtbYB = VNBY;
        
        % Longitudes
        DXA = X1 - XA; DXA3 = X2 - XA;
        DYA = Y1 - YA; DYA3 = Y2 - YA;
        DXB = XB - X2; DXB3 = XB - X1;
        DYB = YB - Y2; DYB3 = YB - Y1;
        
        gt1DXA = gt1X1;  gt2DXA3 = gt2X2;
        gtaDXA = -gtaXA; gtaDXA3 = -gtaXA;
        gt1DYA = gt1Y1;  gt2DYA3 = gt2Y2;
        gtaDYA = -gtaYA; gtaDYA3 = -gtaYA;
        gt2DXB = -gt2X2; gt1DXB3 = -gt1X1;
        gtbDXB = gtbXB;  gtbDXB3 = gtbXB;
        gt2DYB = -gt2Y2; gt1DYB3 = -gt1Y1;
        gtbDYB = gtbYB;  gtbDYB3 = gtbYB;
        
        LQA = DXA*DXA + DYA*DYA;
        LQA3 = DXA3*DXA3 + DYA3*DYA3;
        LQB = DXB*DXB + DYB*DYB;
        LQB3 = DXB3*DXB3 + DYB3*DYB3;
        
        gt1LQA = CD2*DXA*gt1DXA + CD2*DYA*gt1DYA;
        gtaLQA = CD2*DXA*gtaDXA + CD2*DYA*gtaDYA;
        gt2LQA3 = CD2*DXA3*gt2DXA3 + CD2*DYA3*gt2DYA3;
        gtaLQA3 = CD2*DXA3*gtaDXA3 + CD2*DYA3*gtaDYA3;
        gt2LQB = CD2*DXB*gt2DXB + CD2*DYB*gt2DYB;
        gtbLQB = CD2*DXB*gtbDXB + CD2*DYB*gtbDYB;
        gt1LQB3 = CD2*DXB3*gt1DXB3 + CD2*DYB3*gt1DYB3;
        gtbLQB3 = CD2*DXB3*gtbDXB3 + CD2*DYB3*gtbDYB3;
        
        LEA = sqrt(LQA);
        LIA = CD1/LEA;
        LEA3 = sqrt(LQA3);
        LIA3 = CD1/LEA3;
        LEB = sqrt(LQB);
        LIB = CD1/LEB;
        LEB3 = sqrt(LQB3);
        LIB3 = CD1/LEB3;
        
        gt1LEA = CD5*LIA*gt1LQA;
        gtaLEA = CD5*LIA*gtaLQA;
        gt1LIA = -LIA*LIA*gt1LEA;
        gtaLIA = -LIA*LIA*gtaLEA;
        gt2LEA3 = CD5*LIA3*gt2LQA3;
        gtaLEA3 = CD5*LIA3*gtaLQA3;
        gt2LIA3 = -LIA3*LIA3*gt2LEA3;
        gtaLIA3 = -LIA3*LIA3*gtaLEA3;
        gt2LEB = CD5*LIB*gt2LQB;
        gtbLEB = CD5*LIB*gtbLQB;
        gt2LIB = -LIB*LIB*gt2LEB;
        gtbLIB = -LIB*LIB*gtbLEB;
        gt1LEB3 = CD5*LIB3*gt1LQB3;
        gtbLEB3 = CD5*LIB3*gtbLQB3;
        gt1LIB3 = -LIB3*LIB3*gt1LEB3;
        gtbLIB3 = -LIB3*LIB3*gtbLEB3;
        
        SGA = sign(DXA*DY - DYA*DX);
        SGB = sign(DX*DYB - DY*DXB);
        
        LA1 = LEA + LE + LEA3;   LA2 = LEA + LE - LEA3;
        LA3 = LEA3 + LEA - LE;   LA4 = LE + LEA3 - LEA;
        
        gt1LA1 = gt1LEA + gt1LE;    gt1LA2 = gt1LEA + gt1LE;
        gt2LA1 = gt2LE + gt2LEA3;   gt2LA2 = gt2LE - gt2LEA3;
        gtaLA1 = gtaLEA + gtaLEA3;  gtaLA2 = gtaLEA - gtaLEA3;
        gt1LA3 = gt1LEA - gt1LE;    gt1LA4 = gt1LE - gt1LEA;
        gt2LA3 = gt2LEA3 - gt2LE;   gt2LA4 = gt2LE + gt2LEA3;
        gtaLA3 = gtaLEA3 + gtaLEA;  gtaLA4 = gtaLEA3 - gtaLEA;
        
        LA12 = LA1*LA2;          LA34 = LA3*LA4;
        LA14 = LA12*LA34;        LIAS = LI*LIA;
        LAR4 = sqrt(LA14);       LISA = LIAS*LIA3;
        
        gt1LA12 = LA1*gt1LA2 + LA2*gt1LA1;      gt1LA34 = LA3*gt1LA4 + LA4*gt1LA3;
        gt2LA12 = LA1*gt2LA2 + LA2*gt2LA1;      gt2LA34 = LA3*gt2LA4 + LA4*gt2LA3;
        gtaLA12 = LA1*gtaLA2 + LA2*gtaLA1;      gtaLA34 = LA3*gtaLA4 + LA4*gtaLA3;
        gt1LA14 = LA12*gt1LA34 + LA34*gt1LA12;  gt1LIAS = LI*gt1LIA + LIA*gt1LI;
        gt2LA14 = LA12*gt2LA34 + LA34*gt2LA12;  gt2LIAS = LIA*gt2LI;
        gtaLA14 = LA12*gtaLA34 + LA34*gtaLA12;  gtaLIAS = LI*gtaLIA;
        gt1LAR4 = CD5/LAR4*gt1LA14;             gt1LISA = LIA3*gt1LIAS;
        gt2LAR4 = CD5/LAR4*gt2LA14;             gt2LISA = LIAS*gt2LIA3 + LIA3*gt2LIAS;
        gtaLAR4 = CD5/LAR4*gtaLA14;             gtaLISA = LIAS*gtaLIA3 + LIA3*gtaLIAS;
        
        C1 = LAR4*LISA; % Curvatura nodo 1
        C1 = SGA*C1;    % Curvatura nodo 1
        
        gt1C1 = LAR4*gt1LISA + LISA*gt1LAR4;
        gt2C1 = LAR4*gt2LISA + LISA*gt2LAR4;
        gtaC1 = LAR4*gtaLISA + LISA*gtaLAR4;
        gt1C1 = SGA*gt1C1;
        gt2C1 = SGA*gt2C1;
        gtaC1 = SGA*gtaC1;
        
        LB1 = LEB + LE + LEB3;   LB2 = LEB + LE - LEB3;
        LB3 = LEB3 + LEB - LE;   LB4 = LE + LEB3 - LEB;
        
        gt1LB1 = gt1LE + gt1LEB3;   gt1LB2 = gt1LE - gt1LEB3;
        gt2LB1 = gt2LEB + gt2LE;    gt2LB2 = gt2LEB + gt2LE;
        gtbLB1 = gtbLEB + gtbLEB3;  gtbLB2 = gtbLEB - gtbLEB3;
        gt1LB3 = gt1LEB3 - gt1LE;   gt1LB4 = gt1LE + gt1LEB3;
        gt2LB3 = gt2LEB - gt2LE;    gt2LB4 = gt2LE - gt2LEB;
        gtbLB3 = gtbLEB3 + gtbLEB;  gtbLB4 = gtbLEB3 - gtbLEB;
        
        LB12 = LB1*LB2;          LB34 = LB3*LB4;
        LB14 = LB12*LB34;        LIBS = LI*LIB;
        LBR4 = sqrt(LB14);       LISB = LIBS*LIB3;
        
        gt1LB12 = LB1*gt1LB2 + LB2*gt1LB1;      gt1LB34 = LB3*gt1LB4 + LB4*gt1LB3;
        gt2LB12 = LB1*gt2LB2 + LB2*gt2LB1;      gt2LB34 = LB3*gt2LB4 + LB4*gt2LB3;
        gtbLB12 = LB1*gtbLB2 + LB2*gtbLB1;      gtbLB34 = LB3*gtbLB4 + LB4*gtbLB3;
        gt1LB14 = LB12*gt1LB34 + LB34*gt1LB12;  gt1LIBS = LIB*gt1LI;
        gt2LB14 = LB12*gt2LB34 + LB34*gt2LB12;  gt2LIBS = LI*gt2LIB + LIB*gt2LI;
        gtbLB14 = LB12*gtbLB34 + LB34*gtbLB12;  gtbLIBS = LI*gtbLIB;
        gt1LBR4 = CD5/LBR4*gt1LB14;             gt1LISB = LIBS*gt1LIB3 + LIB3*gt1LIBS;
        gt2LBR4 = CD5/LBR4*gt2LB14;             gt2LISB = LIB3*gt2LIBS;
        gtbLBR4 = CD5/LBR4*gtbLB14;             gtbLISB = LIBS*gtbLIB3 + LIB3*gtbLIBS;
        
        C2 = LBR4*LISB; % Curvatura nodo 2
        C2 = SGB*C2;    % Curvatura nodo 2
        
        gt1C2 = LBR4*gt1LISB + LISB*gt1LBR4;
        gt2C2 = LBR4*gt2LISB + LISB*gt2LBR4;
        gtbC2 = LBR4*gtbLISB + LISB*gtbLBR4;
        gt1C2 = SGB*gt1C2;
        gt2C2 = SGB*gt2C2;
        gtbC2 = SGB*gtbC2;
        
        XP1 = XP1 - EST.Gamma*C1;
        XP2 = XP2 - EST.Gamma*C2;
        
        gt1XP1 = gt1XP1 - EST.Gamma*gt1C1;
        gt2XP1 = gt2XP1 - EST.Gamma*gt2C1;
        gtaXP1 = gtaXP1 - EST.Gamma*gtaC1;
        gt1XP2 = gt1XP2 - EST.Gamma*gt1C2;
        gt2XP2 = gt2XP2 - EST.Gamma*gt2C2;
        gtbXP2 = gtbXP2 - EST.Gamma*gtbC2;
        
        NX = +LI*DY;
        NY = -LI*DX;
        
        gt1NX = +LI*gt1DY + DY*gt1LI;
        gt2NX = +LI*gt2DY + DY*gt2LI;
        gt1NY = -LI*gt1DX - DX*gt1LI;
        gt2NY = -LI*gt2DX - DX*gt2LI;
        
        XP1X = XP1*NX;
        XP1Y = XP1*NY;
        XP2X = XP2*NX;
        XP2Y = XP2*NY;
        
        gt1XP1X = XP1*gt1NX + NX*gt1XP1;
        gt2XP1X = XP1*gt2NX + NX*gt2XP1;
        gtaXP1X = NX*gtaXP1;
        gt1XP1Y = XP1*gt1NY + NY*gt1XP1;
        gt2XP1Y = XP1*gt2NY + NY*gt2XP1;
        gtaXP1Y = NY*gtaXP1;
        gt1XP2X = XP2*gt1NX + NX*gt1XP2;
        gt2XP2X = XP2*gt2NX + NX*gt2XP2;
        gtbXP2X = NX*gtbXP2;
        gt1XP2Y = XP2*gt1NY + NY*gt1XP2;
        gt2XP2Y = XP2*gt2NY + NY*gt2XP2;
        gtbXP2Y = NY*gtbXP2;
        % -----------------------------------------------------------------
        
        R1X = X1 - XF;
        R1Y = Y1 - YF;
        R2X = X2 - XF;
        R2Y = Y2 - YF;
        
        gt1R1X = gt1X1;
        gt1R1Y = gt1Y1;
        gt2R2X = gt2X2;
        gt2R2Y = gt2Y2;
        
        Ge11 = CD0;  Ge12 = CD0;
        Ge13 = CD0;  Ge14 = CD0;
        Ge21 = CD0;  Ge22 = CD0;
        Ge23 = CD0;  Ge24 = CD0;
        
        He11 = CD0;  He12 = CD0;
        He13 = CD0;  He14 = CD0;
        He21 = CD0;  He22 = CD0;
        He23 = CD0;  He24 = CD0;
        
        gt1Ge11 = CD0; gt2Ge11 = CD0;
        gt1Ge12 = CD0; gt2Ge12 = CD0;
        gt1Ge13 = CD0; gt2Ge13 = CD0;
        gt1Ge14 = CD0; gt2Ge14 = CD0;
        gt1Ge21 = CD0; gt2Ge21 = CD0;
        gt1Ge22 = CD0; gt2Ge22 = CD0;
        gt1Ge23 = CD0; gt2Ge23 = CD0;
        gt1Ge24 = CD0; gt2Ge24 = CD0;
        
        gt1He11 = CD0; gt2He11 = CD0;
        gt1He12 = CD0; gt2He12 = CD0;
        gt1He13 = CD0; gt2He13 = CD0;
        gt1He14 = CD0; gt2He14 = CD0;
        gt1He21 = CD0; gt2He21 = CD0;
        gt1He22 = CD0; gt2He22 = CD0;
        gt1He23 = CD0; gt2He23 = CD0;
        gt1He24 = CD0; gt2He24 = CD0;
        
        Ra1 = sqrt(R1X*R1X + R1Y*R1Y);
        Ra2 = sqrt(R2X*R2X + R2Y*R2Y);
        
        if (Ra1 < TolRa)
            
            C4L = CM4*LI;
            
            gt1C4L = CM4*gt1LI;
            gt2C4L = CM4*gt2LI;
            
            C4DXX = C4L*DXX;
            C4DXY = C4L*DXY;
            C4DYY = C4L*DYY;
            
            gt1C4DXX = C4L*gt1DXX + DXX*gt1C4L;
            gt2C4DXX = C4L*gt2DXX + DXX*gt2C4L;
            gt1C4DXY = C4L*gt1DXY + DXY*gt1C4L;
            gt2C4DXY = C4L*gt2DXY + DXY*gt2C4L;
            gt1C4DYY = C4L*gt1DYY + DYY*gt1C4L;
            gt2C4DYY = C4L*gt2DYY + DYY*gt2C4L;
            
            CO1 = LogLE - CD6;
            CO2 = LogLE - CD5;
            
            gt1CO1 = gt1LogLE;
            gt2CO1 = gt2LogLE;
            gt1CO2 = gt1LogLE;
            gt2CO2 = gt2LogLE;
            
            C6LC1 = CM6*LE*CO1;
            C6LC2 = CM6*LE*CO2;
            
            gt1C6LC1 = CM6*(LE*gt1CO1 + CO1*gt1LE);
            gt2C6LC1 = CM6*(LE*gt2CO1 + CO1*gt2LE);
            gt1C6LC2 = CM6*(LE*gt1CO2 + CO2*gt1LE);
            gt2C6LC2 = CM6*(LE*gt2CO2 + CO2*gt2LE);
            
            Ge11 = C6LC1 - C4DXX;
            Ge22 = C6LC1 - C4DYY;
            Ge13 = C6LC2 - C4DXX;
            Ge24 = C6LC2 - C4DYY;
            Ge12 = -C4DXY;
            Ge14 = Ge12;
            Ge21 = Ge12;
            Ge23 = Ge12;
            
            gt1Ge11 = gt1C6LC1 - gt1C4DXX;
            gt2Ge11 = gt2C6LC1 - gt2C4DXX;
            gt1Ge22 = gt1C6LC1 - gt1C4DYY;
            gt2Ge22 = gt2C6LC1 - gt2C4DYY;
            gt1Ge13 = gt1C6LC2 - gt1C4DXX;
            gt2Ge13 = gt2C6LC2 - gt2C4DXX;
            gt1Ge24 = gt1C6LC2 - gt1C4DYY;
            gt2Ge24 = gt2C6LC2 - gt2C4DYY;
            gt1Ge12 = -gt1C4DXY;
            gt2Ge12 = -gt2C4DXY;
            gt1Ge14 = gt1Ge12;
            gt2Ge14 = gt2Ge12;
            gt1Ge21 = gt1Ge12;
            gt2Ge21 = gt2Ge12;
            gt1Ge23 = gt1Ge12;
            gt2Ge23 = gt2Ge12;
            
            He14 = CM7;
            He23 = -CM7;
            
        elseif (Ra2 < TolRa)
            
            C4L = CM4*LI;
            
            gt1C4L = CM4*gt1LI;
            gt2C4L = CM4*gt2LI;
            
            C4DXX = C4L*DXX;
            C4DXY = C4L*DXY;
            C4DYY = C4L*DYY;
            
            gt1C4DXX = C4L*gt1DXX + DXX*gt1C4L;
            gt2C4DXX = C4L*gt2DXX + DXX*gt2C4L;
            gt1C4DXY = C4L*gt1DXY + DXY*gt1C4L;
            gt2C4DXY = C4L*gt2DXY + DXY*gt2C4L;
            gt1C4DYY = C4L*gt1DYY + DYY*gt1C4L;
            gt2C4DYY = C4L*gt2DYY + DYY*gt2C4L;
            
            CO1 = LogLE - CD5;
            CO2 = LogLE - CD6;
            
            gt1CO1 = gt1LogLE;
            gt2CO1 = gt2LogLE;
            gt1CO2 = gt1LogLE;
            gt2CO2 = gt2LogLE;
            
            C6LC1 = CM6*LE*CO1;
            C6LC2 = CM6*LE*CO2;
            
            gt1C6LC1 = CM6*(LE*gt1CO1 + CO1*gt1LE);
            gt2C6LC1 = CM6*(LE*gt2CO1 + CO1*gt2LE);
            gt1C6LC2 = CM6*(LE*gt1CO2 + CO2*gt1LE);
            gt2C6LC2 = CM6*(LE*gt2CO2 + CO2*gt2LE);
            
            Ge11 = C6LC1 - C4DXX;
            Ge22 = C6LC1 - C4DYY;
            Ge13 = C6LC2 - C4DXX;
            Ge24 = C6LC2 - C4DYY;
            Ge12 = -C4DXY;
            Ge14 = Ge12;
            Ge21 = Ge12;
            Ge23 = Ge12;
            
            gt1Ge11 = gt1C6LC1 - gt1C4DXX;
            gt2Ge11 = gt2C6LC1 - gt2C4DXX;
            gt1Ge22 = gt1C6LC1 - gt1C4DYY;
            gt2Ge22 = gt2C6LC1 - gt2C4DYY;
            gt1Ge13 = gt1C6LC2 - gt1C4DXX;
            gt2Ge13 = gt2C6LC2 - gt2C4DXX;
            gt1Ge24 = gt1C6LC2 - gt1C4DYY;
            gt2Ge24 = gt2C6LC2 - gt2C4DYY;
            gt1Ge12 = -gt1C4DXY;
            gt2Ge12 = -gt2C4DXY;
            gt1Ge14 = gt1Ge12;
            gt2Ge14 = gt2Ge12;
            gt1Ge21 = gt1Ge12;
            gt2Ge21 = gt2Ge12;
            gt1Ge23 = gt1Ge12;
            gt2Ge23 = gt2Ge12;
            
            He12 = -CM7;
            He21 = CM7;
            
        else
            
            for PI = 1:NG
                
                N1 = CD5*(CD1 - XI(PI));
                N2 = CD5*(CD1 + XI(PI));
                
                RX = N1*R1X + N2*R2X;
                RY = N1*R1Y + N2*R2Y;
                
                gt1RX = N1*gt1R1X;
                gt2RX = N2*gt2R2X;
                gt1RY = N1*gt1R1Y;
                gt2RY = N2*gt2R2Y;
                
                RQ = RX*RX + RY*RY;
                RE = sqrt(RQ);
                RI = CD1/RE;
                LogRE = log(RE);
                
                gt1RQ = CD2*(RX*gt1RX + RY*gt1RY);
                gt2RQ = CD2*(RX*gt2RX + RY*gt2RY);
                gt1RE = CD5*RI*gt1RQ;
                gt2RE = CD5*RI*gt2RQ;
                gt1RI = -RI*RI*gt1RE;
                gt2RI = -RI*RI*gt2RE;
                gt1LogRE = RI*gt1RE;
                gt2LogRE = RI*gt2RE;
                
                dRX = RI*RX;
                dRY = RI*RY;
                
                gt1dRX = RI*gt1RX + RX*gt1RI;
                gt2dRX = RI*gt2RX + RX*gt2RI;
                gt1dRY = RI*gt1RY + RY*gt1RI;
                gt2dRY = RI*gt2RY + RY*gt2RI;
                
                dRN = dRX*DY - dRY*DX;
                dRT = dRX*DX + dRY*DY;
                
                gt1dRN = dRX*gt1DY + DY*gt1dRX - dRY*gt1DX - DX*gt1dRY;
                gt2dRN = dRX*gt2DY + DY*gt2dRX - dRY*gt2DX - DX*gt2dRY;
                gt1dRT = dRX*gt1DX + DX*gt1dRX + dRY*gt1DY + DY*gt1dRY;
                gt2dRT = dRX*gt2DX + DX*gt2dRX + dRY*gt2DY + DY*gt2dRY;
                
                C6R = CM5*LogRE;
                
                gt1C6R = CM5*gt1LogRE;
                gt2C6R = CM5*gt2LogRE;
                
                dRXX = dRX*dRX;
                dRXY = dRX*dRY;
                dRYY = dRY*dRY;
                
                gt1dRXX = CD2*dRX*gt1dRX;
                gt2dRXX = CD2*dRX*gt2dRX;
                gt1dRXY = dRX*gt1dRY + dRY*gt1dRX;
                gt2dRXY = dRX*gt2dRY + dRY*gt2dRX;
                gt1dRYY = CD2*dRY*gt1dRY;
                gt2dRYY = CD2*dRY*gt2dRY;
                
                UXX = C6R - CM3*dRXX;
                UYY = C6R - CM3*dRYY;
                UXY = -CM3*dRXY;
                
                gt1UXX = gt1C6R - CM3*gt1dRXX;
                gt2UXX = gt2C6R - CM3*gt2dRXX;
                gt1UYY = gt1C6R - CM3*gt1dRYY;
                gt2UYY = gt2C6R - CM3*gt2dRYY;
                gt1UXY = -CM3*gt1dRXY;
                gt2UXY = -CM3*gt2dRXY;
                
                C4XX = CM7 + CM2*dRXX;
                C4YY = CM7 + CM2*dRYY;
                C4XY = CM2*dRXY;
                
                gt1C4XX = CM2*gt1dRXX;
                gt2C4XX = CM2*gt2dRXX;
                gt1C4YY = CM2*gt1dRYY;
                gt2C4YY = CM2*gt2dRYY;
                gt1C4XY = CM2*gt1dRXY;
                gt2C4XY = CM2*gt2dRXY;
                
                dRNRI = RI*dRN;
                dRTRI = RI*dRT;
                
                gt1dRNRI = RI*gt1dRN + dRN*gt1RI;
                gt2dRNRI = RI*gt2dRN + dRN*gt2RI;
                gt1dRTRI = RI*gt1dRT + dRT*gt1RI;
                gt2dRTRI = RI*gt2dRT + dRT*gt2RI;
                
                PXX = C4XX*dRNRI;
                PYY = C4YY*dRNRI;
                
                gt1PXX = C4XX*gt1dRNRI + dRNRI*gt1C4XX;
                gt2PXX = C4XX*gt2dRNRI + dRNRI*gt2C4XX;
                gt1PYY = C4YY*gt1dRNRI + dRNRI*gt1C4YY;
                gt2PYY = C4YY*gt2dRNRI + dRNRI*gt2C4YY;
                
                PXY = C4XY*dRNRI + CM7*dRTRI;
                PYX = C4XY*dRNRI - CM7*dRTRI;
                
                gt1PXY = C4XY*gt1dRNRI + dRNRI*gt1C4XY + CM7*gt1dRTRI;
                gt2PXY = C4XY*gt2dRNRI + dRNRI*gt2C4XY + CM7*gt2dRTRI;
                gt1PYX = C4XY*gt1dRNRI + dRNRI*gt1C4XY - CM7*gt1dRTRI;
                gt2PYX = C4XY*gt2dRNRI + dRNRI*gt2C4XY - CM7*gt2dRTRI;
                
                W5 = CD5*WI(PI);
                
                WUXX = W5*LE*UXX;
                WUXY = W5*LE*UXY;
                WUYY = W5*LE*UYY;
                
                gt1WUXX = W5*(LE*gt1UXX + UXX*gt1LE);
                gt2WUXX = W5*(LE*gt2UXX + UXX*gt2LE);
                gt1WUXY = W5*(LE*gt1UXY + UXY*gt1LE);
                gt2WUXY = W5*(LE*gt2UXY + UXY*gt2LE);
                gt1WUYY = W5*(LE*gt1UYY + UYY*gt1LE);
                gt2WUYY = W5*(LE*gt2UYY + UYY*gt2LE);
                
                WPXX = W5*PXX;
                WPXY = W5*PXY;
                WPYX = W5*PYX;
                WPYY = W5*PYY;
                
                gt1WPXX = W5*gt1PXX;
                gt2WPXX = W5*gt2PXX;
                gt1WPXY = W5*gt1PXY;
                gt2WPXY = W5*gt2PXY;
                gt1WPYX = W5*gt1PYX;
                gt2WPYX = W5*gt2PYX;
                gt1WPYY = W5*gt1PYY;
                gt2WPYY = W5*gt2PYY;
                
                He11 = He11 + WPXX*N1;
                He12 = He12 + WPXY*N1;
                He13 = He13 + WPXX*N2;
                He14 = He14 + WPXY*N2;
                
                He21 = He21 + WPYX*N1;
                He22 = He22 + WPYY*N1;
                He23 = He23 + WPYX*N2;
                He24 = He24 + WPYY*N2;
                
                gt1He11 = gt1He11 + gt1WPXX*N1;
                gt2He11 = gt2He11 + gt2WPXX*N1;
                gt1He12 = gt1He12 + gt1WPXY*N1;
                gt2He12 = gt2He12 + gt2WPXY*N1;
                gt1He13 = gt1He13 + gt1WPXX*N2;
                gt2He13 = gt2He13 + gt2WPXX*N2;
                gt1He14 = gt1He14 + gt1WPXY*N2;
                gt2He14 = gt2He14 + gt2WPXY*N2;
                
                gt1He21 = gt1He21 + gt1WPYX*N1;
                gt2He21 = gt2He21 + gt2WPYX*N1;
                gt1He22 = gt1He22 + gt1WPYY*N1;
                gt2He22 = gt2He22 + gt2WPYY*N1;
                gt1He23 = gt1He23 + gt1WPYX*N2;
                gt2He23 = gt2He23 + gt2WPYX*N2;
                gt1He24 = gt1He24 + gt1WPYY*N2;
                gt2He24 = gt2He24 + gt2WPYY*N2;
                
                Ge11 = Ge11 + WUXX*N1;
                Ge12 = Ge12 + WUXY*N1;
                Ge13 = Ge13 + WUXX*N2;
                Ge14 = Ge14 + WUXY*N2;
                
                Ge21 = Ge21 + WUXY*N1;
                Ge22 = Ge22 + WUYY*N1;
                Ge23 = Ge23 + WUXY*N2;
                Ge24 = Ge24 + WUYY*N2;
                
                gt1Ge11 = gt1Ge11 + gt1WUXX*N1;
                gt2Ge11 = gt2Ge11 + gt2WUXX*N1;
                gt1Ge12 = gt1Ge12 + gt1WUXY*N1;
                gt2Ge12 = gt2Ge12 + gt2WUXY*N1;
                gt1Ge13 = gt1Ge13 + gt1WUXX*N2;
                gt2Ge13 = gt2Ge13 + gt2WUXX*N2;
                gt1Ge14 = gt1Ge14 + gt1WUXY*N2;
                gt2Ge14 = gt2Ge14 + gt2WUXY*N2;
                
                gt1Ge21 = gt1Ge21 + gt1WUXY*N1;
                gt2Ge21 = gt2Ge21 + gt2WUXY*N1;
                gt1Ge22 = gt1Ge22 + gt1WUYY*N1;
                gt2Ge22 = gt2Ge22 + gt2WUYY*N1;
                gt1Ge23 = gt1Ge23 + gt1WUXY*N2;
                gt2Ge23 = gt2Ge23 + gt2WUXY*N2;
                gt1Ge24 = gt1Ge24 + gt1WUYY*N2;
                gt2Ge24 = gt2Ge24 + gt2WUYY*N2;
                
            end
        end
        
        % Monto matrices MatH y MatG --------------------------------------
        
        % MatH(GrF,GrU) = MatH(GrF,GrU) + He;
        % MatG(GrF,GrP) = MatG(GrF,GrP) + Ge;
        
        % G(ResF1) = G(ResF1) - Ge11*XP1X - Ge12*XP1Y - Ge13*XP2X - Ge14*XP2Y;
        % G(ResF2) = G(ResF2) - Ge21*XP1X - Ge22*XP1Y - Ge23*XP2X - Ge24*XP2Y;
        % G(ResF1) = G(ResF1) + He11*XU1X + He12*XU1Y + He13*XU2X + He14*XU2Y;
        % G(ResF2) = G(ResF2) + He21*XU1X + He22*XU1Y + He23*XU2X + He24*XU2Y;
        
        gG(ResF1,VarT1) = gG(ResF1,VarT1) - XP1X*gt1Ge11 - XP1Y*gt1Ge12 - XP2X*gt1Ge13 - XP2Y*gt1Ge14;
        gG(ResF1,VarT1) = gG(ResF1,VarT1) - Ge11*gt1XP1X - Ge12*gt1XP1Y - Ge13*gt1XP2X - Ge14*gt1XP2Y;
        gG(ResF1,VarT2) = gG(ResF1,VarT2) - XP1X*gt2Ge11 - XP1Y*gt2Ge12 - XP2X*gt2Ge13 - XP2Y*gt2Ge14;
        gG(ResF1,VarT2) = gG(ResF1,VarT2) - Ge11*gt2XP1X - Ge12*gt2XP1Y - Ge13*gt2XP2X - Ge14*gt2XP2Y;
        gG(ResF1,VarTA) = gG(ResF1,VarTA) - Ge11*gtaXP1X - Ge12*gtaXP1Y;
        gG(ResF1,VarTB) = gG(ResF1,VarTB) - Ge13*gtbXP2X - Ge14*gtbXP2Y;
        gG(ResF2,VarT1) = gG(ResF2,VarT1) - XP1X*gt1Ge21 - XP1Y*gt1Ge22 - XP2X*gt1Ge23 - XP2Y*gt1Ge24;
        gG(ResF2,VarT1) = gG(ResF2,VarT1) - Ge21*gt1XP1X - Ge22*gt1XP1Y - Ge23*gt1XP2X - Ge24*gt1XP2Y;
        gG(ResF2,VarT2) = gG(ResF2,VarT2) - XP1X*gt2Ge21 - XP1Y*gt2Ge22 - XP2X*gt2Ge23 - XP2Y*gt2Ge24;
        gG(ResF2,VarT2) = gG(ResF2,VarT2) - Ge21*gt2XP1X - Ge22*gt2XP1Y - Ge23*gt2XP2X - Ge24*gt2XP2Y;
        gG(ResF2,VarTA) = gG(ResF2,VarTA) - Ge21*gtaXP1X - Ge22*gtaXP1Y;
        gG(ResF2,VarTB) = gG(ResF2,VarTB) - Ge23*gtbXP2X - Ge24*gtbXP2Y;
        gG(ResF1,VarT1) = gG(ResF1,VarT1) + XU1X*gt1He11 + XU1Y*gt1He12 + XU2X*gt1He13 + XU2Y*gt1He14;
        gG(ResF1,VarT2) = gG(ResF1,VarT2) + XU1X*gt2He11 + XU1Y*gt2He12 + XU2X*gt2He13 + XU2Y*gt2He14;
        gG(ResF2,VarT1) = gG(ResF2,VarT1) + XU1X*gt1He21 + XU1Y*gt1He22 + XU2X*gt1He23 + XU2Y*gt1He24;
        gG(ResF2,VarT2) = gG(ResF2,VarT2) + XU1X*gt2He21 + XU1Y*gt2He22 + XU2X*gt2He23 + XU2Y*gt2He24;
        
        gG(ResF1,VarU1X) = gG(ResF1,VarU1X) + He11*gu1xXU1X;
        gG(ResF1,VarU1Y) = gG(ResF1,VarU1Y) + He12*gu1yXU1Y;
        gG(ResF1,VarU2X) = gG(ResF1,VarU2X) + He13*gu2xXU2X;
        gG(ResF1,VarU2Y) = gG(ResF1,VarU2Y) + He14*gu2yXU2Y;
        gG(ResF2,VarU1X) = gG(ResF2,VarU1X) + He21*gu1xXU1X;
        gG(ResF2,VarU1Y) = gG(ResF2,VarU1Y) + He22*gu1yXU1Y;
        gG(ResF2,VarU2X) = gG(ResF2,VarU2X) + He23*gu2xXU2X;
        gG(ResF2,VarU2Y) = gG(ResF2,VarU2Y) + He24*gu2yXU2Y;
        
        % Desplazamientos rígidos -----------------------------------------
        
        % MatH(GrF,GrF) = MatH(GrF,GrF) - He(1:2,1:2) - He(1:2,3:4);
        
        % G(ResF1) = G(ResF1) - He11*XUFX - He12*XUFY - He13*XUFX - He14*XUFY;
        % G(ResF2) = G(ResF2) - He21*XUFX - He22*XUFY - He23*XUFX - He24*XUFY;
        
        gG(ResF1,VarT1) = gG(ResF1,VarT1) - XUFX*gt1He11 - XUFY*gt1He12 - XUFX*gt1He13 - XUFY*gt1He14;
        gG(ResF1,VarT2) = gG(ResF1,VarT2) - XUFX*gt2He11 - XUFY*gt2He12 - XUFX*gt2He13 - XUFY*gt2He14;
        gG(ResF2,VarT1) = gG(ResF2,VarT1) - XUFX*gt1He21 - XUFY*gt1He22 - XUFX*gt1He23 - XUFY*gt1He24;
        gG(ResF2,VarT2) = gG(ResF2,VarT2) - XUFX*gt2He21 - XUFY*gt2He22 - XUFX*gt2He23 - XUFY*gt2He24;
        
    end
    
    % Variable lamda en el sistema ----------------------------------------
    VarLAM = MatXL;
    
    % XLAM = X(VarLAM);
    
    % G(ResF1) = G(ResF1) + EST.PN(2*NodF-1)*XLAM;
    % G(ResF2) = G(ResF2) + EST.PN(2*NodF)*XLAM;
    
    gG(ResF1,VarLAM) = gG(ResF1,VarLAM) + EST.PN(2*NodF-1);
    gG(ResF2,VarLAM) = gG(ResF2,VarLAM) + EST.PN(2*NodF);

    VarPFX = MatXP(2*NodF-1);
    VarPFY = MatXP(2*NodF);
    
    % XPFX = X(VarPFX);
    % XPFY = X(VarPFY);
    
    % G(ResLAM) = G(ResLAM) + EST.PN(2*NodF-1)*XPFX;
    % G(ResLAM) = G(ResLAM) + EST.PN(2*NodF)*XPFY;

    gG(ResLAM,VarPFX) = gG(ResLAM,VarPFX) + EST.PN(2*NodF-1);
    gG(ResLAM,VarPFY) = gG(ResLAM,VarPFY) + EST.PN(2*NodF);
    
end

for NodF = 1:EST.NumNod2
    
    VarTF = MatXT(NodF);
    
    VNFX = EST.VZ(NodF,1);
    VNFY = EST.VZ(NodF,2);
    
    XF = EST.MatNod2(NodF,1) + VNFX*X(VarTF);
    YF = EST.MatNod2(NodF,2) + VNFY*X(VarTF);
    
    gtfXF = VNFX;
    gtfYF = VNFY;
    
    VarUFX = MatXU(2*NodF-1);
    VarUFY = MatXU(2*NodF);
    
    XUFX = X(VarUFX);
    XUFY = X(VarUFY);
    
    gufxXUFX = CD1;
    gufyXUFY = CD1;
    
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
        
        gp1xXP1X = CD1;
        gp1yXP1Y = CD1;
        gp2xXP2X = CD1;
        gp2yXP2Y = CD1;
        
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
        
        gtfR1X = -gtfXF;
        gtfR1Y = -gtfYF;
        gtfR2X = -gtfXF;
        gtfR2Y = -gtfYF;
        
        Ge11 = CD0;  Ge12 = CD0;
        Ge13 = CD0;  Ge14 = CD0;
        Ge21 = CD0;  Ge22 = CD0;
        Ge23 = CD0;  Ge24 = CD0;
        
        He11 = CD0;  He12 = CD0;
        He13 = CD0;  He14 = CD0;
        He21 = CD0;  He22 = CD0;
        He23 = CD0;  He24 = CD0;
        
        gtfGe11 = CD0;
        gtfGe12 = CD0;
        gtfGe13 = CD0;
        gtfGe14 = CD0;
        gtfGe21 = CD0;
        gtfGe22 = CD0;
        gtfGe23 = CD0;
        gtfGe24 = CD0;
        
        gtfHe11 = CD0;
        gtfHe12 = CD0;
        gtfHe13 = CD0;
        gtfHe14 = CD0;
        gtfHe21 = CD0;
        gtfHe22 = CD0;
        gtfHe23 = CD0;
        gtfHe24 = CD0;
        
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
                
                gtfRX = N1*gtfR1X + N2*gtfR2X;
                gtfRY = N1*gtfR1Y + N2*gtfR2Y;
                
                RQ = RX*RX + RY*RY;
                RE = sqrt(RQ);
                RI = CD1/RE;
                LogRE = log(RE);
                
                gtfRQ = CD2*(RX*gtfRX + RY*gtfRY);
                gtfRE = CD5*RI*gtfRQ;
                gtfRI = -RI*RI*gtfRE;
                gtfLogRE = RI*gtfRE;
                
                dRX = RI*RX;
                dRY = RI*RY;
                
                gtfdRX = RI*gtfRX + RX*gtfRI;
                gtfdRY = RI*gtfRY + RY*gtfRI;
                
                dRN = dRX*DY - dRY*DX;
                dRT = dRX*DX + dRY*DY;
                
                gtfdRN = DY*gtfdRX - DX*gtfdRY;
                gtfdRT = DX*gtfdRX + DY*gtfdRY;
                
                C6R = CM5*LogRE;
                
                gtfC6R = CM5*gtfLogRE;
                
                dRXX = dRX*dRX;
                dRXY = dRX*dRY;
                dRYY = dRY*dRY;
                
                gtfdRXX = CD2*dRX*gtfdRX;
                gtfdRXY = dRX*gtfdRY + dRY*gtfdRX;
                gtfdRYY = CD2*dRY*gtfdRY;
                
                UXX = C6R - CM3*dRXX;
                UYY = C6R - CM3*dRYY;
                UXY = -CM3*dRXY;
                
                gtfUXX = gtfC6R - CM3*gtfdRXX;
                gtfUYY = gtfC6R - CM3*gtfdRYY;
                gtfUXY = -CM3*gtfdRXY;
                
                C4XX = CM7 + CM2*dRXX;
                C4YY = CM7 + CM2*dRYY;
                C4XY = CM2*dRXY;
                
                gtfC4XX = CM2*gtfdRXX;
                gtfC4YY = CM2*gtfdRYY;
                gtfC4XY = CM2*gtfdRXY;
                
                dRNRI = RI*dRN;
                dRTRI = RI*dRT;
                
                gtfdRNRI = RI*gtfdRN + dRN*gtfRI;
                gtfdRTRI = RI*gtfdRT + dRT*gtfRI;
                
                PXX = C4XX*dRNRI;
                PYY = C4YY*dRNRI;
                
                gtfPXX = C4XX*gtfdRNRI + dRNRI*gtfC4XX;
                gtfPYY = C4YY*gtfdRNRI + dRNRI*gtfC4YY;
                
                PXY = C4XY*dRNRI + CM7*dRTRI;
                PYX = C4XY*dRNRI - CM7*dRTRI;
                
                gtfPXY = C4XY*gtfdRNRI + dRNRI*gtfC4XY + CM7*gtfdRTRI;
                gtfPYX = C4XY*gtfdRNRI + dRNRI*gtfC4XY - CM7*gtfdRTRI;
                
                W5 = CD5*WI(PI);
                
                WUXX = W5*LE*UXX;
                WUXY = W5*LE*UXY;
                WUYY = W5*LE*UYY;
                
                gtfWUXX = W5*LE*gtfUXX;
                gtfWUXY = W5*LE*gtfUXY;
                gtfWUYY = W5*LE*gtfUYY;
                
                WPXX = W5*PXX;
                WPXY = W5*PXY;
                WPYX = W5*PYX;
                WPYY = W5*PYY;
                
                gtfWPXX = W5*gtfPXX;
                gtfWPXY = W5*gtfPXY;
                gtfWPYX = W5*gtfPYX;
                gtfWPYY = W5*gtfPYY;
                
                He11 = He11 + WPXX*N1;
                He12 = He12 + WPXY*N1;
                He13 = He13 + WPXX*N2;
                He14 = He14 + WPXY*N2;
                
                He21 = He21 + WPYX*N1;
                He22 = He22 + WPYY*N1;
                He23 = He23 + WPYX*N2;
                He24 = He24 + WPYY*N2;
                
                gtfHe11 = gtfHe11 + gtfWPXX*N1;
                gtfHe12 = gtfHe12 + gtfWPXY*N1;
                gtfHe13 = gtfHe13 + gtfWPXX*N2;
                gtfHe14 = gtfHe14 + gtfWPXY*N2;
                
                gtfHe21 = gtfHe21 + gtfWPYX*N1;
                gtfHe22 = gtfHe22 + gtfWPYY*N1;
                gtfHe23 = gtfHe23 + gtfWPYX*N2;
                gtfHe24 = gtfHe24 + gtfWPYY*N2;
                
                Ge11 = Ge11 + WUXX*N1;
                Ge12 = Ge12 + WUXY*N1;
                Ge13 = Ge13 + WUXX*N2;
                Ge14 = Ge14 + WUXY*N2;
                
                Ge21 = Ge21 + WUXY*N1;
                Ge22 = Ge22 + WUYY*N1;
                Ge23 = Ge23 + WUXY*N2;
                Ge24 = Ge24 + WUYY*N2;
                
                gtfGe11 = gtfGe11 + gtfWUXX*N1;
                gtfGe12 = gtfGe12 + gtfWUXY*N1;
                gtfGe13 = gtfGe13 + gtfWUXX*N2;
                gtfGe14 = gtfGe14 + gtfWUXY*N2;
                
                gtfGe21 = gtfGe21 + gtfWUXY*N1;
                gtfGe22 = gtfGe22 + gtfWUYY*N1;
                gtfGe23 = gtfGe23 + gtfWUXY*N2;
                gtfGe24 = gtfGe24 + gtfWUYY*N2;
                
            end
        end
        
        % Monto matrices MatH y MatG --------------------------------------
        
        % MatH(GrF,GrU) = MatH(GrF,GrU) + He;
        % MatG(GrF,GrP) = MatG(GrF,GrP) + Ge;
        
        % G(ResF1) = G(ResF1) - Ge11*XP1X - Ge12*XP1Y - Ge13*XP2X - Ge14*XP2Y;
        % G(ResF2) = G(ResF2) - Ge21*XP1X - Ge22*XP1Y - Ge23*XP2X - Ge24*XP2Y;
        % G(ResF1) = G(ResF1) + He11*XU1X + He12*XU1Y + He13*XU2X + He14*XU2Y;
        % G(ResF2) = G(ResF2) + He21*XU1X + He22*XU1Y + He23*XU2X + He24*XU2Y;
        
        gG(ResF1,VarTF) = gG(ResF1,VarTF) - XP1X*gtfGe11 - XP1Y*gtfGe12 - XP2X*gtfGe13 - XP2Y*gtfGe14;
        gG(ResF2,VarTF) = gG(ResF2,VarTF) - XP1X*gtfGe21 - XP1Y*gtfGe22 - XP2X*gtfGe23 - XP2Y*gtfGe24;
        gG(ResF1,VarTF) = gG(ResF1,VarTF) + XU1X*gtfHe11 + XU1Y*gtfHe12 + XU2X*gtfHe13 + XU2Y*gtfHe14;
        gG(ResF2,VarTF) = gG(ResF2,VarTF) + XU1X*gtfHe21 + XU1Y*gtfHe22 + XU2X*gtfHe23 + XU2Y*gtfHe24;
        
        gG(ResF1,VarP1X) = gG(ResF1,VarP1X) - Ge11*gp1xXP1X;
        gG(ResF1,VarP1Y) = gG(ResF1,VarP1Y) - Ge12*gp1yXP1Y;
        gG(ResF1,VarP2X) = gG(ResF1,VarP2X) - Ge13*gp2xXP2X;
        gG(ResF1,VarP2Y) = gG(ResF1,VarP2Y) - Ge14*gp2yXP2Y;
        gG(ResF2,VarP1X) = gG(ResF2,VarP1X) - Ge21*gp1xXP1X;
        gG(ResF2,VarP1Y) = gG(ResF2,VarP1Y) - Ge22*gp1yXP1Y;
        gG(ResF2,VarP2X) = gG(ResF2,VarP2X) - Ge23*gp2xXP2X;
        gG(ResF2,VarP2Y) = gG(ResF2,VarP2Y) - Ge24*gp2yXP2Y;
        
        % Desplazamientos rígidos -----------------------------------------
        
        % MatH(GrF,GrF) = MatH(GrF,GrF) - He(1:2,1:2) - He(1:2,3:4);
        
        % G(ResF1) = G(ResF1) - He11*XUFX - He12*XUFY - He13*XUFX - He14*XUFY;
        % G(ResF2) = G(ResF2) - He21*XUFX - He22*XUFY - He23*XUFX - He24*XUFY;
        
        gG(ResF1,VarTF) = gG(ResF1,VarTF) - XUFX*gtfHe11 - XUFY*gtfHe12 - XUFX*gtfHe13 - XUFY*gtfHe14;
        gG(ResF2,VarTF) = gG(ResF2,VarTF) - XUFX*gtfHe21 - XUFY*gtfHe22 - XUFX*gtfHe23 - XUFY*gtfHe24;
        
        gG(ResF1,VarUFX) = gG(ResF1,VarUFX) - He11*gufxXUFX - He13*gufxXUFX;
        gG(ResF1,VarUFY) = gG(ResF1,VarUFY) - He12*gufyXUFY - He14*gufyXUFY;
        gG(ResF2,VarUFX) = gG(ResF2,VarUFX) - He21*gufxXUFX - He23*gufxXUFX;
        gG(ResF2,VarUFY) = gG(ResF2,VarUFY) - He22*gufyXUFY - He24*gufyXUFY;
        
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
        
        gu1xXU1X = CD1;
        gu1yXU1Y = CD1;
        gu2xXU2X = CD1;
        gu2yXU2Y = CD1;
        
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
        
        gt1X1 = VN1X;
        gt1Y1 = VN1Y;
        gt2X2 = VN2X;
        gt2Y2 = VN2Y;
        
        DX = X2 - X1;
        DY = Y2 - Y1;
        
        gt1DX = -gt1X1;
        gt2DX =  gt2X2;
        gt1DY = -gt1Y1;
        gt2DY =  gt2Y2;
        
        DXX = DX*DX;
        DXY = DX*DY;
        DYY = DY*DY;
        
        gt1DXX = CD2*DX*gt1DX;
        gt2DXX = CD2*DX*gt2DX;
        gt1DXY = DX*gt1DY + DY*gt1DX;
        gt2DXY = DX*gt2DY + DY*gt2DX;
        gt1DYY = CD2*DY*gt1DY;
        gt2DYY = CD2*DY*gt2DY;
        
        LQ = DXX + DYY;
        LE = sqrt(LQ);
        LI = CD1/LE;
        LogLE = log(LE);
        
        gt1LQ = gt1DXX + gt1DYY;
        gt2LQ = gt2DXX + gt2DYY;
        gt1LE = CD5*LI*gt1LQ;
        gt2LE = CD5*LI*gt2LQ;
        gt1LI = -LI*LI*gt1LE;
        gt2LI = -LI*LI*gt2LE;
        gt1LogLE = LI*gt1LE;
        gt2LogLE = LI*gt2LE;
        
        % Presion exterior ------------------------------------------------
        XP1 = -EST.p0; % presion impuesta
        XP2 = -EST.p0; % presion impuesta
        
        gt1XP1 = CD0;
        gt2XP1 = CD0;
        gtaXP1 = CD0;
        gt1XP2 = CD0;
        gt2XP2 = CD0;
        gtbXP2 = CD0;
        
        XP1 = XP1 + EST.bA*Y1; % resto solucion particular
        XP2 = XP2 + EST.bA*Y2; % resto solucion particular
        
        gt1XP1 = gt1XP1 + EST.bA*gt1Y1;
        gt2XP2 = gt2XP2 + EST.bA*gt2Y2;
        
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
        
        gtaXA = VNAX;
        gtaYA = VNAY;
        gtbXB = VNBX;
        gtbYB = VNBY;
        
        % Longitudes
        DXA = X1 - XA; DXA3 = X2 - XA;
        DYA = Y1 - YA; DYA3 = Y2 - YA;
        DXB = XB - X2; DXB3 = XB - X1;
        DYB = YB - Y2; DYB3 = YB - Y1;
        
        gt1DXA = gt1X1;  gt2DXA3 = gt2X2;
        gtaDXA = -gtaXA; gtaDXA3 = -gtaXA;
        gt1DYA = gt1Y1;  gt2DYA3 = gt2Y2;
        gtaDYA = -gtaYA; gtaDYA3 = -gtaYA;
        gt2DXB = -gt2X2; gt1DXB3 = -gt1X1;
        gtbDXB = gtbXB;  gtbDXB3 = gtbXB;
        gt2DYB = -gt2Y2; gt1DYB3 = -gt1Y1;
        gtbDYB = gtbYB;  gtbDYB3 = gtbYB;
        
        LQA = DXA*DXA + DYA*DYA;
        LQA3 = DXA3*DXA3 + DYA3*DYA3;
        LQB = DXB*DXB + DYB*DYB;
        LQB3 = DXB3*DXB3 + DYB3*DYB3;
        
        gt1LQA = CD2*DXA*gt1DXA + CD2*DYA*gt1DYA;
        gtaLQA = CD2*DXA*gtaDXA + CD2*DYA*gtaDYA;
        gt2LQA3 = CD2*DXA3*gt2DXA3 + CD2*DYA3*gt2DYA3;
        gtaLQA3 = CD2*DXA3*gtaDXA3 + CD2*DYA3*gtaDYA3;
        gt2LQB = CD2*DXB*gt2DXB + CD2*DYB*gt2DYB;
        gtbLQB = CD2*DXB*gtbDXB + CD2*DYB*gtbDYB;
        gt1LQB3 = CD2*DXB3*gt1DXB3 + CD2*DYB3*gt1DYB3;
        gtbLQB3 = CD2*DXB3*gtbDXB3 + CD2*DYB3*gtbDYB3;
        
        LEA = sqrt(LQA);
        LIA = CD1/LEA;
        LEA3 = sqrt(LQA3);
        LIA3 = CD1/LEA3;
        LEB = sqrt(LQB);
        LIB = CD1/LEB;
        LEB3 = sqrt(LQB3);
        LIB3 = CD1/LEB3;
        
        gt1LEA = CD5*LIA*gt1LQA;
        gtaLEA = CD5*LIA*gtaLQA;
        gt1LIA = -LIA*LIA*gt1LEA;
        gtaLIA = -LIA*LIA*gtaLEA;
        gt2LEA3 = CD5*LIA3*gt2LQA3;
        gtaLEA3 = CD5*LIA3*gtaLQA3;
        gt2LIA3 = -LIA3*LIA3*gt2LEA3;
        gtaLIA3 = -LIA3*LIA3*gtaLEA3;
        gt2LEB = CD5*LIB*gt2LQB;
        gtbLEB = CD5*LIB*gtbLQB;
        gt2LIB = -LIB*LIB*gt2LEB;
        gtbLIB = -LIB*LIB*gtbLEB;
        gt1LEB3 = CD5*LIB3*gt1LQB3;
        gtbLEB3 = CD5*LIB3*gtbLQB3;
        gt1LIB3 = -LIB3*LIB3*gt1LEB3;
        gtbLIB3 = -LIB3*LIB3*gtbLEB3;
        
        SGA = sign(DXA*DY - DYA*DX);
        SGB = sign(DX*DYB - DY*DXB);
        
        LA1 = LEA + LE + LEA3;   LA2 = LEA + LE - LEA3;
        LA3 = LEA3 + LEA - LE;   LA4 = LE + LEA3 - LEA;
        
        gt1LA1 = gt1LEA + gt1LE;    gt1LA2 = gt1LEA + gt1LE;
        gt2LA1 = gt2LE + gt2LEA3;   gt2LA2 = gt2LE - gt2LEA3;
        gtaLA1 = gtaLEA + gtaLEA3;  gtaLA2 = gtaLEA - gtaLEA3;
        gt1LA3 = gt1LEA - gt1LE;    gt1LA4 = gt1LE - gt1LEA;
        gt2LA3 = gt2LEA3 - gt2LE;   gt2LA4 = gt2LE + gt2LEA3;
        gtaLA3 = gtaLEA3 + gtaLEA;  gtaLA4 = gtaLEA3 - gtaLEA;
        
        LA12 = LA1*LA2;          LA34 = LA3*LA4;
        LA14 = LA12*LA34;        LIAS = LI*LIA;
        LAR4 = sqrt(LA14);       LISA = LIAS*LIA3;
        
        gt1LA12 = LA1*gt1LA2 + LA2*gt1LA1;      gt1LA34 = LA3*gt1LA4 + LA4*gt1LA3;
        gt2LA12 = LA1*gt2LA2 + LA2*gt2LA1;      gt2LA34 = LA3*gt2LA4 + LA4*gt2LA3;
        gtaLA12 = LA1*gtaLA2 + LA2*gtaLA1;      gtaLA34 = LA3*gtaLA4 + LA4*gtaLA3;
        gt1LA14 = LA12*gt1LA34 + LA34*gt1LA12;  gt1LIAS = LI*gt1LIA + LIA*gt1LI;
        gt2LA14 = LA12*gt2LA34 + LA34*gt2LA12;  gt2LIAS = LIA*gt2LI;
        gtaLA14 = LA12*gtaLA34 + LA34*gtaLA12;  gtaLIAS = LI*gtaLIA;
        gt1LAR4 = CD5/LAR4*gt1LA14;             gt1LISA = LIA3*gt1LIAS;
        gt2LAR4 = CD5/LAR4*gt2LA14;             gt2LISA = LIAS*gt2LIA3 + LIA3*gt2LIAS;
        gtaLAR4 = CD5/LAR4*gtaLA14;             gtaLISA = LIAS*gtaLIA3 + LIA3*gtaLIAS;
        
        C1 = LAR4*LISA; % Curvatura nodo 1
        C1 = SGA*C1;    % Curvatura nodo 1
        
        gt1C1 = LAR4*gt1LISA + LISA*gt1LAR4;
        gt2C1 = LAR4*gt2LISA + LISA*gt2LAR4;
        gtaC1 = LAR4*gtaLISA + LISA*gtaLAR4;
        gt1C1 = SGA*gt1C1;
        gt2C1 = SGA*gt2C1;
        gtaC1 = SGA*gtaC1;
        
        LB1 = LEB + LE + LEB3;   LB2 = LEB + LE - LEB3;
        LB3 = LEB3 + LEB - LE;   LB4 = LE + LEB3 - LEB;
        
        gt1LB1 = gt1LE + gt1LEB3;   gt1LB2 = gt1LE - gt1LEB3;
        gt2LB1 = gt2LEB + gt2LE;    gt2LB2 = gt2LEB + gt2LE;
        gtbLB1 = gtbLEB + gtbLEB3;  gtbLB2 = gtbLEB - gtbLEB3;
        gt1LB3 = gt1LEB3 - gt1LE;   gt1LB4 = gt1LE + gt1LEB3;
        gt2LB3 = gt2LEB - gt2LE;    gt2LB4 = gt2LE - gt2LEB;
        gtbLB3 = gtbLEB3 + gtbLEB;  gtbLB4 = gtbLEB3 - gtbLEB;
        
        LB12 = LB1*LB2;          LB34 = LB3*LB4;
        LB14 = LB12*LB34;        LIBS = LI*LIB;
        LBR4 = sqrt(LB14);       LISB = LIBS*LIB3;
        
        gt1LB12 = LB1*gt1LB2 + LB2*gt1LB1;      gt1LB34 = LB3*gt1LB4 + LB4*gt1LB3;
        gt2LB12 = LB1*gt2LB2 + LB2*gt2LB1;      gt2LB34 = LB3*gt2LB4 + LB4*gt2LB3;
        gtbLB12 = LB1*gtbLB2 + LB2*gtbLB1;      gtbLB34 = LB3*gtbLB4 + LB4*gtbLB3;
        gt1LB14 = LB12*gt1LB34 + LB34*gt1LB12;  gt1LIBS = LIB*gt1LI;
        gt2LB14 = LB12*gt2LB34 + LB34*gt2LB12;  gt2LIBS = LI*gt2LIB + LIB*gt2LI;
        gtbLB14 = LB12*gtbLB34 + LB34*gtbLB12;  gtbLIBS = LI*gtbLIB;
        gt1LBR4 = CD5/LBR4*gt1LB14;             gt1LISB = LIBS*gt1LIB3 + LIB3*gt1LIBS;
        gt2LBR4 = CD5/LBR4*gt2LB14;             gt2LISB = LIB3*gt2LIBS;
        gtbLBR4 = CD5/LBR4*gtbLB14;             gtbLISB = LIBS*gtbLIB3 + LIB3*gtbLIBS;
        
        C2 = LBR4*LISB; % Curvatura nodo 2
        C2 = SGB*C2;    % Curvatura nodo 2
        
        gt1C2 = LBR4*gt1LISB + LISB*gt1LBR4;
        gt2C2 = LBR4*gt2LISB + LISB*gt2LBR4;
        gtbC2 = LBR4*gtbLISB + LISB*gtbLBR4;
        gt1C2 = SGB*gt1C2;
        gt2C2 = SGB*gt2C2;
        gtbC2 = SGB*gtbC2;
        
        XP1 = XP1 - EST.Gamma*C1;
        XP2 = XP2 - EST.Gamma*C2;
        
        gt1XP1 = gt1XP1 - EST.Gamma*gt1C1;
        gt2XP1 = gt2XP1 - EST.Gamma*gt2C1;
        gtaXP1 = gtaXP1 - EST.Gamma*gtaC1;
        gt1XP2 = gt1XP2 - EST.Gamma*gt1C2;
        gt2XP2 = gt2XP2 - EST.Gamma*gt2C2;
        gtbXP2 = gtbXP2 - EST.Gamma*gtbC2;
        
        NX = +LI*DY;
        NY = -LI*DX;
        
        gt1NX = +LI*gt1DY + DY*gt1LI;
        gt2NX = +LI*gt2DY + DY*gt2LI;
        gt1NY = -LI*gt1DX - DX*gt1LI;
        gt2NY = -LI*gt2DX - DX*gt2LI;
        
        XP1X = XP1*NX;
        XP1Y = XP1*NY;
        XP2X = XP2*NX;
        XP2Y = XP2*NY;
        
        gt1XP1X = XP1*gt1NX + NX*gt1XP1;
        gt2XP1X = XP1*gt2NX + NX*gt2XP1;
        gtaXP1X = NX*gtaXP1;
        gt1XP1Y = XP1*gt1NY + NY*gt1XP1;
        gt2XP1Y = XP1*gt2NY + NY*gt2XP1;
        gtaXP1Y = NY*gtaXP1;
        gt1XP2X = XP2*gt1NX + NX*gt1XP2;
        gt2XP2X = XP2*gt2NX + NX*gt2XP2;
        gtbXP2X = NX*gtbXP2;
        gt1XP2Y = XP2*gt1NY + NY*gt1XP2;
        gt2XP2Y = XP2*gt2NY + NY*gt2XP2;
        gtbXP2Y = NY*gtbXP2;
        % -----------------------------------------------------------------
        
        R1X = X1 - XF;
        R1Y = Y1 - YF;
        R2X = X2 - XF;
        R2Y = Y2 - YF;
        
        gt1R1X = gt1X1;
        gtfR1X = -gtfXF;
        gt1R1Y = gt1Y1;
        gtfR1Y = -gtfYF;
        gt2R2X = gt2X2;
        gtfR2X = -gtfXF;
        gt2R2Y = gt2Y2;
        gtfR2Y = -gtfYF;
        
        Ge11 = CD0;  Ge12 = CD0;
        Ge13 = CD0;  Ge14 = CD0;
        Ge21 = CD0;  Ge22 = CD0;
        Ge23 = CD0;  Ge24 = CD0;
        
        He11 = CD0;  He12 = CD0;
        He13 = CD0;  He14 = CD0;
        He21 = CD0;  He22 = CD0;
        He23 = CD0;  He24 = CD0;
        
        gt1Ge11 = CD0; gt2Ge11 = CD0; gtfGe11 = CD0;
        gt1Ge12 = CD0; gt2Ge12 = CD0; gtfGe12 = CD0;
        gt1Ge13 = CD0; gt2Ge13 = CD0; gtfGe13 = CD0;
        gt1Ge14 = CD0; gt2Ge14 = CD0; gtfGe14 = CD0;
        gt1Ge21 = CD0; gt2Ge21 = CD0; gtfGe21 = CD0;
        gt1Ge22 = CD0; gt2Ge22 = CD0; gtfGe22 = CD0;
        gt1Ge23 = CD0; gt2Ge23 = CD0; gtfGe23 = CD0;
        gt1Ge24 = CD0; gt2Ge24 = CD0; gtfGe24 = CD0;
        
        gt1He11 = CD0; gt2He11 = CD0; gtfHe11 = CD0;
        gt1He12 = CD0; gt2He12 = CD0; gtfHe12 = CD0;
        gt1He13 = CD0; gt2He13 = CD0; gtfHe13 = CD0;
        gt1He14 = CD0; gt2He14 = CD0; gtfHe14 = CD0;
        gt1He21 = CD0; gt2He21 = CD0; gtfHe21 = CD0;
        gt1He22 = CD0; gt2He22 = CD0; gtfHe22 = CD0;
        gt1He23 = CD0; gt2He23 = CD0; gtfHe23 = CD0;
        gt1He24 = CD0; gt2He24 = CD0; gtfHe24 = CD0;
        
        Ra1 = sqrt(R1X*R1X + R1Y*R1Y);
        Ra2 = sqrt(R2X*R2X + R2Y*R2Y);
        
        if (Ra1 < TolRa)
            
            C4L = CM4*LI;
            
            gt1C4L = CM4*gt1LI;
            gt2C4L = CM4*gt2LI;
            
            C4DXX = C4L*DXX;
            C4DXY = C4L*DXY;
            C4DYY = C4L*DYY;
            
            gt1C4DXX = C4L*gt1DXX + DXX*gt1C4L;
            gt2C4DXX = C4L*gt2DXX + DXX*gt2C4L;
            gt1C4DXY = C4L*gt1DXY + DXY*gt1C4L;
            gt2C4DXY = C4L*gt2DXY + DXY*gt2C4L;
            gt1C4DYY = C4L*gt1DYY + DYY*gt1C4L;
            gt2C4DYY = C4L*gt2DYY + DYY*gt2C4L;
            
            CO1 = LogLE - CD6;
            CO2 = LogLE - CD5;
            
            gt1CO1 = gt1LogLE;
            gt2CO1 = gt2LogLE;
            gt1CO2 = gt1LogLE;
            gt2CO2 = gt2LogLE;
            
            C6LC1 = CM6*LE*CO1;
            C6LC2 = CM6*LE*CO2;
            
            gt1C6LC1 = CM6*(LE*gt1CO1 + CO1*gt1LE);
            gt2C6LC1 = CM6*(LE*gt2CO1 + CO1*gt2LE);
            gt1C6LC2 = CM6*(LE*gt1CO2 + CO2*gt1LE);
            gt2C6LC2 = CM6*(LE*gt2CO2 + CO2*gt2LE);
            
            Ge11 = C6LC1 - C4DXX;
            Ge22 = C6LC1 - C4DYY;
            Ge13 = C6LC2 - C4DXX;
            Ge24 = C6LC2 - C4DYY;
            Ge12 = -C4DXY;
            Ge14 = Ge12;
            Ge21 = Ge12;
            Ge23 = Ge12;
            
            gt1Ge11 = gt1C6LC1 - gt1C4DXX;
            gt2Ge11 = gt2C6LC1 - gt2C4DXX;
            gt1Ge22 = gt1C6LC1 - gt1C4DYY;
            gt2Ge22 = gt2C6LC1 - gt2C4DYY;
            gt1Ge13 = gt1C6LC2 - gt1C4DXX;
            gt2Ge13 = gt2C6LC2 - gt2C4DXX;
            gt1Ge24 = gt1C6LC2 - gt1C4DYY;
            gt2Ge24 = gt2C6LC2 - gt2C4DYY;
            gt1Ge12 = -gt1C4DXY;
            gt2Ge12 = -gt2C4DXY;
            gt1Ge14 = gt1Ge12;
            gt2Ge14 = gt2Ge12;
            gt1Ge21 = gt1Ge12;
            gt2Ge21 = gt2Ge12;
            gt1Ge23 = gt1Ge12;
            gt2Ge23 = gt2Ge12;
            
            He14 = CM7;
            He23 = -CM7;
            
        elseif (Ra2 < TolRa)
            
            C4L = CM4*LI;
            
            gt1C4L = CM4*gt1LI;
            gt2C4L = CM4*gt2LI;
            
            C4DXX = C4L*DXX;
            C4DXY = C4L*DXY;
            C4DYY = C4L*DYY;
            
            gt1C4DXX = C4L*gt1DXX + DXX*gt1C4L;
            gt2C4DXX = C4L*gt2DXX + DXX*gt2C4L;
            gt1C4DXY = C4L*gt1DXY + DXY*gt1C4L;
            gt2C4DXY = C4L*gt2DXY + DXY*gt2C4L;
            gt1C4DYY = C4L*gt1DYY + DYY*gt1C4L;
            gt2C4DYY = C4L*gt2DYY + DYY*gt2C4L;
            
            CO1 = LogLE - CD5;
            CO2 = LogLE - CD6;
            
            gt1CO1 = gt1LogLE;
            gt2CO1 = gt2LogLE;
            gt1CO2 = gt1LogLE;
            gt2CO2 = gt2LogLE;
            
            C6LC1 = CM6*LE*CO1;
            C6LC2 = CM6*LE*CO2;
            
            gt1C6LC1 = CM6*(LE*gt1CO1 + CO1*gt1LE);
            gt2C6LC1 = CM6*(LE*gt2CO1 + CO1*gt2LE);
            gt1C6LC2 = CM6*(LE*gt1CO2 + CO2*gt1LE);
            gt2C6LC2 = CM6*(LE*gt2CO2 + CO2*gt2LE);
            
            Ge11 = C6LC1 - C4DXX;
            Ge22 = C6LC1 - C4DYY;
            Ge13 = C6LC2 - C4DXX;
            Ge24 = C6LC2 - C4DYY;
            Ge12 = -C4DXY;
            Ge14 = Ge12;
            Ge21 = Ge12;
            Ge23 = Ge12;
            
            gt1Ge11 = gt1C6LC1 - gt1C4DXX;
            gt2Ge11 = gt2C6LC1 - gt2C4DXX;
            gt1Ge22 = gt1C6LC1 - gt1C4DYY;
            gt2Ge22 = gt2C6LC1 - gt2C4DYY;
            gt1Ge13 = gt1C6LC2 - gt1C4DXX;
            gt2Ge13 = gt2C6LC2 - gt2C4DXX;
            gt1Ge24 = gt1C6LC2 - gt1C4DYY;
            gt2Ge24 = gt2C6LC2 - gt2C4DYY;
            gt1Ge12 = -gt1C4DXY;
            gt2Ge12 = -gt2C4DXY;
            gt1Ge14 = gt1Ge12;
            gt2Ge14 = gt2Ge12;
            gt1Ge21 = gt1Ge12;
            gt2Ge21 = gt2Ge12;
            gt1Ge23 = gt1Ge12;
            gt2Ge23 = gt2Ge12;
            
            He12 = -CM7;
            He21 = CM7;
            
        else
            
            for PI = 1:NG
                
                N1 = CD5*(CD1 - XI(PI));
                N2 = CD5*(CD1 + XI(PI));
                
                RX = N1*R1X + N2*R2X;
                RY = N1*R1Y + N2*R2Y;
                
                gt1RX = N1*gt1R1X;
                gt2RX = N2*gt2R2X;
                gtfRX = N1*gtfR1X + N2*gtfR2X;
                gt1RY = N1*gt1R1Y;
                gt2RY = N2*gt2R2Y;
                gtfRY = N1*gtfR1Y + N2*gtfR2Y;
                
                RQ = RX*RX + RY*RY;
                RE = sqrt(RQ);
                RI = CD1/RE;
                LogRE = log(RE);
                
                gt1RQ = CD2*(RX*gt1RX + RY*gt1RY);
                gt2RQ = CD2*(RX*gt2RX + RY*gt2RY);
                gtfRQ = CD2*(RX*gtfRX + RY*gtfRY);
                gt1RE = CD5*RI*gt1RQ;
                gt2RE = CD5*RI*gt2RQ;
                gtfRE = CD5*RI*gtfRQ;
                gt1RI = -RI*RI*gt1RE;
                gt2RI = -RI*RI*gt2RE;
                gtfRI = -RI*RI*gtfRE;
                gt1LogRE = RI*gt1RE;
                gt2LogRE = RI*gt2RE;
                gtfLogRE = RI*gtfRE;
                
                dRX = RI*RX;
                dRY = RI*RY;
                
                gt1dRX = RI*gt1RX + RX*gt1RI;
                gt2dRX = RI*gt2RX + RX*gt2RI;
                gtfdRX = RI*gtfRX + RX*gtfRI;
                gt1dRY = RI*gt1RY + RY*gt1RI;
                gt2dRY = RI*gt2RY + RY*gt2RI;
                gtfdRY = RI*gtfRY + RY*gtfRI;
                
                dRN = dRX*DY - dRY*DX;
                dRT = dRX*DX + dRY*DY;
                
                gt1dRN = dRX*gt1DY + DY*gt1dRX - dRY*gt1DX - DX*gt1dRY;
                gt2dRN = dRX*gt2DY + DY*gt2dRX - dRY*gt2DX - DX*gt2dRY;
                gtfdRN = DY*gtfdRX - DX*gtfdRY;
                gt1dRT = dRX*gt1DX + DX*gt1dRX + dRY*gt1DY + DY*gt1dRY;
                gt2dRT = dRX*gt2DX + DX*gt2dRX + dRY*gt2DY + DY*gt2dRY;
                gtfdRT = DX*gtfdRX + DY*gtfdRY;
                
                C6R = CM5*LogRE;
                
                gt1C6R = CM5*gt1LogRE;
                gt2C6R = CM5*gt2LogRE;
                gtfC6R = CM5*gtfLogRE;
                
                dRXX = dRX*dRX;
                dRXY = dRX*dRY;
                dRYY = dRY*dRY;
                
                gt1dRXX = CD2*dRX*gt1dRX;
                gt2dRXX = CD2*dRX*gt2dRX;
                gtfdRXX = CD2*dRX*gtfdRX;
                gt1dRXY = dRX*gt1dRY + dRY*gt1dRX;
                gt2dRXY = dRX*gt2dRY + dRY*gt2dRX;
                gtfdRXY = dRX*gtfdRY + dRY*gtfdRX;
                gt1dRYY = CD2*dRY*gt1dRY;
                gt2dRYY = CD2*dRY*gt2dRY;
                gtfdRYY = CD2*dRY*gtfdRY;
                
                UXX = C6R - CM3*dRXX;
                UYY = C6R - CM3*dRYY;
                UXY = -CM3*dRXY;
                
                gt1UXX = gt1C6R - CM3*gt1dRXX;
                gt2UXX = gt2C6R - CM3*gt2dRXX;
                gtfUXX = gtfC6R - CM3*gtfdRXX;
                gt1UYY = gt1C6R - CM3*gt1dRYY;
                gt2UYY = gt2C6R - CM3*gt2dRYY;
                gtfUYY = gtfC6R - CM3*gtfdRYY;
                gt1UXY = -CM3*gt1dRXY;
                gt2UXY = -CM3*gt2dRXY;
                gtfUXY = -CM3*gtfdRXY;
                
                C4XX = CM7 + CM2*dRXX;
                C4YY = CM7 + CM2*dRYY;
                C4XY = CM2*dRXY;
                
                gt1C4XX = CM2*gt1dRXX;
                gt2C4XX = CM2*gt2dRXX;
                gtfC4XX = CM2*gtfdRXX;
                gt1C4YY = CM2*gt1dRYY;
                gt2C4YY = CM2*gt2dRYY;
                gtfC4YY = CM2*gtfdRYY;
                gt1C4XY = CM2*gt1dRXY;
                gt2C4XY = CM2*gt2dRXY;
                gtfC4XY = CM2*gtfdRXY;
                
                dRNRI = RI*dRN;
                dRTRI = RI*dRT;
                
                gt1dRNRI = RI*gt1dRN + dRN*gt1RI;
                gt2dRNRI = RI*gt2dRN + dRN*gt2RI;
                gtfdRNRI = RI*gtfdRN + dRN*gtfRI;
                gt1dRTRI = RI*gt1dRT + dRT*gt1RI;
                gt2dRTRI = RI*gt2dRT + dRT*gt2RI;
                gtfdRTRI = RI*gtfdRT + dRT*gtfRI;
                
                PXX = C4XX*dRNRI;
                PYY = C4YY*dRNRI;
                
                gt1PXX = C4XX*gt1dRNRI + dRNRI*gt1C4XX;
                gt2PXX = C4XX*gt2dRNRI + dRNRI*gt2C4XX;
                gtfPXX = C4XX*gtfdRNRI + dRNRI*gtfC4XX;
                gt1PYY = C4YY*gt1dRNRI + dRNRI*gt1C4YY;
                gt2PYY = C4YY*gt2dRNRI + dRNRI*gt2C4YY;
                gtfPYY = C4YY*gtfdRNRI + dRNRI*gtfC4YY;
                
                PXY = C4XY*dRNRI + CM7*dRTRI;
                PYX = C4XY*dRNRI - CM7*dRTRI;
                
                gt1PXY = C4XY*gt1dRNRI + dRNRI*gt1C4XY + CM7*gt1dRTRI;
                gt2PXY = C4XY*gt2dRNRI + dRNRI*gt2C4XY + CM7*gt2dRTRI;
                gtfPXY = C4XY*gtfdRNRI + dRNRI*gtfC4XY + CM7*gtfdRTRI;
                gt1PYX = C4XY*gt1dRNRI + dRNRI*gt1C4XY - CM7*gt1dRTRI;
                gt2PYX = C4XY*gt2dRNRI + dRNRI*gt2C4XY - CM7*gt2dRTRI;
                gtfPYX = C4XY*gtfdRNRI + dRNRI*gtfC4XY - CM7*gtfdRTRI;
                
                W5 = CD5*WI(PI);
                
                WUXX = W5*LE*UXX;
                WUXY = W5*LE*UXY;
                WUYY = W5*LE*UYY;
                
                gt1WUXX = W5*(LE*gt1UXX + UXX*gt1LE);
                gt2WUXX = W5*(LE*gt2UXX + UXX*gt2LE);
                gtfWUXX = W5*(LE*gtfUXX);
                gt1WUXY = W5*(LE*gt1UXY + UXY*gt1LE);
                gt2WUXY = W5*(LE*gt2UXY + UXY*gt2LE);
                gtfWUXY = W5*(LE*gtfUXY);
                gt1WUYY = W5*(LE*gt1UYY + UYY*gt1LE);
                gt2WUYY = W5*(LE*gt2UYY + UYY*gt2LE);
                gtfWUYY = W5*(LE*gtfUYY);
                
                WPXX = W5*PXX;
                WPXY = W5*PXY;
                WPYX = W5*PYX;
                WPYY = W5*PYY;
                
                gt1WPXX = W5*gt1PXX;
                gt2WPXX = W5*gt2PXX;
                gtfWPXX = W5*gtfPXX;
                gt1WPXY = W5*gt1PXY;
                gt2WPXY = W5*gt2PXY;
                gtfWPXY = W5*gtfPXY;
                gt1WPYX = W5*gt1PYX;
                gt2WPYX = W5*gt2PYX;
                gtfWPYX = W5*gtfPYX;
                gt1WPYY = W5*gt1PYY;
                gt2WPYY = W5*gt2PYY;
                gtfWPYY = W5*gtfPYY;
                
                He11 = He11 + WPXX*N1;
                He12 = He12 + WPXY*N1;
                He13 = He13 + WPXX*N2;
                He14 = He14 + WPXY*N2;
                
                He21 = He21 + WPYX*N1;
                He22 = He22 + WPYY*N1;
                He23 = He23 + WPYX*N2;
                He24 = He24 + WPYY*N2;
                
                gt1He11 = gt1He11 + gt1WPXX*N1;
                gt2He11 = gt2He11 + gt2WPXX*N1;
                gtfHe11 = gtfHe11 + gtfWPXX*N1;
                gt1He12 = gt1He12 + gt1WPXY*N1;
                gt2He12 = gt2He12 + gt2WPXY*N1;
                gtfHe12 = gtfHe12 + gtfWPXY*N1;
                gt1He13 = gt1He13 + gt1WPXX*N2;
                gt2He13 = gt2He13 + gt2WPXX*N2;
                gtfHe13 = gtfHe13 + gtfWPXX*N2;
                gt1He14 = gt1He14 + gt1WPXY*N2;
                gt2He14 = gt2He14 + gt2WPXY*N2;
                gtfHe14 = gtfHe14 + gtfWPXY*N2;
                
                gt1He21 = gt1He21 + gt1WPYX*N1;
                gt2He21 = gt2He21 + gt2WPYX*N1;
                gtfHe21 = gtfHe21 + gtfWPYX*N1;
                gt1He22 = gt1He22 + gt1WPYY*N1;
                gt2He22 = gt2He22 + gt2WPYY*N1;
                gtfHe22 = gtfHe22 + gtfWPYY*N1;
                gt1He23 = gt1He23 + gt1WPYX*N2;
                gt2He23 = gt2He23 + gt2WPYX*N2;
                gtfHe23 = gtfHe23 + gtfWPYX*N2;
                gt1He24 = gt1He24 + gt1WPYY*N2;
                gt2He24 = gt2He24 + gt2WPYY*N2;
                gtfHe24 = gtfHe24 + gtfWPYY*N2;
                
                Ge11 = Ge11 + WUXX*N1;
                Ge12 = Ge12 + WUXY*N1;
                Ge13 = Ge13 + WUXX*N2;
                Ge14 = Ge14 + WUXY*N2;
                
                Ge21 = Ge21 + WUXY*N1;
                Ge22 = Ge22 + WUYY*N1;
                Ge23 = Ge23 + WUXY*N2;
                Ge24 = Ge24 + WUYY*N2;
                
                gt1Ge11 = gt1Ge11 + gt1WUXX*N1;
                gt2Ge11 = gt2Ge11 + gt2WUXX*N1;
                gtfGe11 = gtfGe11 + gtfWUXX*N1;
                gt1Ge12 = gt1Ge12 + gt1WUXY*N1;
                gt2Ge12 = gt2Ge12 + gt2WUXY*N1;
                gtfGe12 = gtfGe12 + gtfWUXY*N1;
                gt1Ge13 = gt1Ge13 + gt1WUXX*N2;
                gt2Ge13 = gt2Ge13 + gt2WUXX*N2;
                gtfGe13 = gtfGe13 + gtfWUXX*N2;
                gt1Ge14 = gt1Ge14 + gt1WUXY*N2;
                gt2Ge14 = gt2Ge14 + gt2WUXY*N2;
                gtfGe14 = gtfGe14 + gtfWUXY*N2;
                
                gt1Ge21 = gt1Ge21 + gt1WUXY*N1;
                gt2Ge21 = gt2Ge21 + gt2WUXY*N1;
                gtfGe21 = gtfGe21 + gtfWUXY*N1;
                gt1Ge22 = gt1Ge22 + gt1WUYY*N1;
                gt2Ge22 = gt2Ge22 + gt2WUYY*N1;
                gtfGe22 = gtfGe22 + gtfWUYY*N1;
                gt1Ge23 = gt1Ge23 + gt1WUXY*N2;
                gt2Ge23 = gt2Ge23 + gt2WUXY*N2;
                gtfGe23 = gtfGe23 + gtfWUXY*N2;
                gt1Ge24 = gt1Ge24 + gt1WUYY*N2;
                gt2Ge24 = gt2Ge24 + gt2WUYY*N2;
                gtfGe24 = gtfGe24 + gtfWUYY*N2;
                
            end
        end
        
        % Monto matrices MatH y MatG --------------------------------------
        
        % MatH(GrF,GrU) = MatH(GrF,GrU) + He;
        % MatG(GrF,GrP) = MatG(GrF,GrP) + Ge;
        
        % G(ResF1) = G(ResF1) - Ge11*XP1X - Ge12*XP1Y - Ge13*XP2X - Ge14*XP2Y;
        % G(ResF2) = G(ResF2) - Ge21*XP1X - Ge22*XP1Y - Ge23*XP2X - Ge24*XP2Y;
        % G(ResF1) = G(ResF1) + He11*XU1X + He12*XU1Y + He13*XU2X + He14*XU2Y;
        % G(ResF2) = G(ResF2) + He21*XU1X + He22*XU1Y + He23*XU2X + He24*XU2Y;
        
        gG(ResF1,VarT1) = gG(ResF1,VarT1) - XP1X*gt1Ge11 - XP1Y*gt1Ge12 - XP2X*gt1Ge13 - XP2Y*gt1Ge14;
        gG(ResF1,VarT1) = gG(ResF1,VarT1) - Ge11*gt1XP1X - Ge12*gt1XP1Y - Ge13*gt1XP2X - Ge14*gt1XP2Y;
        gG(ResF1,VarT2) = gG(ResF1,VarT2) - XP1X*gt2Ge11 - XP1Y*gt2Ge12 - XP2X*gt2Ge13 - XP2Y*gt2Ge14;
        gG(ResF1,VarT2) = gG(ResF1,VarT2) - Ge11*gt2XP1X - Ge12*gt2XP1Y - Ge13*gt2XP2X - Ge14*gt2XP2Y;
        gG(ResF1,VarTA) = gG(ResF1,VarTA) - Ge11*gtaXP1X - Ge12*gtaXP1Y;
        gG(ResF1,VarTB) = gG(ResF1,VarTB) - Ge13*gtbXP2X - Ge14*gtbXP2Y;
        gG(ResF1,VarTF) = gG(ResF1,VarTF) - XP1X*gtfGe11 - XP1Y*gtfGe12 - XP2X*gtfGe13 - XP2Y*gtfGe14;
        gG(ResF2,VarT1) = gG(ResF2,VarT1) - XP1X*gt1Ge21 - XP1Y*gt1Ge22 - XP2X*gt1Ge23 - XP2Y*gt1Ge24;
        gG(ResF2,VarT1) = gG(ResF2,VarT1) - Ge21*gt1XP1X - Ge22*gt1XP1Y - Ge23*gt1XP2X - Ge24*gt1XP2Y;
        gG(ResF2,VarT2) = gG(ResF2,VarT2) - XP1X*gt2Ge21 - XP1Y*gt2Ge22 - XP2X*gt2Ge23 - XP2Y*gt2Ge24;
        gG(ResF2,VarT2) = gG(ResF2,VarT2) - Ge21*gt2XP1X - Ge22*gt2XP1Y - Ge23*gt2XP2X - Ge24*gt2XP2Y;
        gG(ResF2,VarTA) = gG(ResF2,VarTA) - Ge21*gtaXP1X - Ge22*gtaXP1Y;
        gG(ResF2,VarTB) = gG(ResF2,VarTB) - Ge23*gtbXP2X - Ge24*gtbXP2Y;
        gG(ResF2,VarTF) = gG(ResF2,VarTF) - XP1X*gtfGe21 - XP1Y*gtfGe22 - XP2X*gtfGe23 - XP2Y*gtfGe24;
        gG(ResF1,VarT1) = gG(ResF1,VarT1) + XU1X*gt1He11 + XU1Y*gt1He12 + XU2X*gt1He13 + XU2Y*gt1He14;
        gG(ResF1,VarT2) = gG(ResF1,VarT2) + XU1X*gt2He11 + XU1Y*gt2He12 + XU2X*gt2He13 + XU2Y*gt2He14;
        gG(ResF1,VarTF) = gG(ResF1,VarTF) + XU1X*gtfHe11 + XU1Y*gtfHe12 + XU2X*gtfHe13 + XU2Y*gtfHe14;
        gG(ResF2,VarT1) = gG(ResF2,VarT1) + XU1X*gt1He21 + XU1Y*gt1He22 + XU2X*gt1He23 + XU2Y*gt1He24;
        gG(ResF2,VarT2) = gG(ResF2,VarT2) + XU1X*gt2He21 + XU1Y*gt2He22 + XU2X*gt2He23 + XU2Y*gt2He24;
        gG(ResF2,VarTF) = gG(ResF2,VarTF) + XU1X*gtfHe21 + XU1Y*gtfHe22 + XU2X*gtfHe23 + XU2Y*gtfHe24;
        
        gG(ResF1,VarU1X) = gG(ResF1,VarU1X) + He11*gu1xXU1X;
        gG(ResF1,VarU1Y) = gG(ResF1,VarU1Y) + He12*gu1yXU1Y;
        gG(ResF1,VarU2X) = gG(ResF1,VarU2X) + He13*gu2xXU2X;
        gG(ResF1,VarU2Y) = gG(ResF1,VarU2Y) + He14*gu2yXU2Y;
        gG(ResF2,VarU1X) = gG(ResF2,VarU1X) + He21*gu1xXU1X;
        gG(ResF2,VarU1Y) = gG(ResF2,VarU1Y) + He22*gu1yXU1Y;
        gG(ResF2,VarU2X) = gG(ResF2,VarU2X) + He23*gu2xXU2X;
        gG(ResF2,VarU2Y) = gG(ResF2,VarU2Y) + He24*gu2yXU2Y;
        
        % Desplazamientos rígidos -----------------------------------------
        
        % MatH(GrF,GrF) = MatH(GrF,GrF) - He(1:2,1:2) - He(1:2,3:4);
        
        % G(ResF1) = G(ResF1) - He11*XUFX - He12*XUFY - He13*XUFX - He14*XUFY;
        % G(ResF2) = G(ResF2) - He21*XUFX - He22*XUFY - He23*XUFX - He24*XUFY;
        
        gG(ResF1,VarT1) = gG(ResF1,VarT1) - XUFX*gt1He11 - XUFY*gt1He12 - XUFX*gt1He13 - XUFY*gt1He14;
        gG(ResF1,VarT2) = gG(ResF1,VarT2) - XUFX*gt2He11 - XUFY*gt2He12 - XUFX*gt2He13 - XUFY*gt2He14;
        gG(ResF1,VarTF) = gG(ResF1,VarTF) - XUFX*gtfHe11 - XUFY*gtfHe12 - XUFX*gtfHe13 - XUFY*gtfHe14;
        gG(ResF2,VarT1) = gG(ResF2,VarT1) - XUFX*gt1He21 - XUFY*gt1He22 - XUFX*gt1He23 - XUFY*gt1He24;
        gG(ResF2,VarT2) = gG(ResF2,VarT2) - XUFX*gt2He21 - XUFY*gt2He22 - XUFX*gt2He23 - XUFY*gt2He24;
        gG(ResF2,VarTF) = gG(ResF2,VarTF) - XUFX*gtfHe21 - XUFY*gtfHe22 - XUFX*gtfHe23 - XUFY*gtfHe24;
        
        gG(ResF1,VarUFX) = gG(ResF1,VarUFX) - He11*gufxXUFX - He13*gufxXUFX;
        gG(ResF1,VarUFY) = gG(ResF1,VarUFY) - He12*gufyXUFY - He14*gufyXUFY;
        gG(ResF2,VarUFX) = gG(ResF2,VarUFX) - He21*gufxXUFX - He23*gufxXUFX;
        gG(ResF2,VarUFY) = gG(ResF2,VarUFY) - He22*gufyXUFY - He24*gufyXUFY;
        
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
    
    gt1X1 = VN1X;
    gt1Y1 = VN1Y;
    gt2X2 = VN2X;
    gt2Y2 = VN2Y;
    
    VarU1X = MatXU(2*Nod1-1);
    VarU1Y = MatXU(2*Nod1);
    VarU2X = MatXU(2*Nod2-1);
    VarU2Y = MatXU(2*Nod2);
    
    % Solucion particular
    XU1XP = CD0;
    XU1YP = CD0;
    XU2XP = CD0;
    XU2YP = CD0;
    
    gt1XU1XP = CD0;
    gt1XU1YP = CD0;
    gt2XU2XP = CD0;
    gt2XU2YP = CD0;

    XU1X = X(VarU1X) + XU1XP; % sumo solucion particular
    XU1Y = X(VarU1Y) + XU1YP; % sumo solucion particular
    XU2X = X(VarU2X) + XU2XP; % sumo solucion particular
    XU2Y = X(VarU2Y) + XU2YP; % sumo solucion particular
    
    gt1XU1X = gt1XU1XP;
    gu1xXU1X = CD1;
    gt1XU1Y = gt1XU1YP;
    gu1yXU1Y = CD1;
    gt2XU2X = gt2XU2XP;
    gu2xXU2X = CD1;
    gt2XU2Y = gt2XU2YP;
    gu2yXU2Y = CD1;
    
    DX = X2 - X1;
    DY = Y2 - Y1;
    
    gt1DX = -gt1X1;
    gt2DX =  gt2X2;
    gt1DY = -gt1Y1;
    gt2DY =  gt2Y2;
    
    % U1N = XU1X*DY - XU1Y*DX;
    % U2N = XU2X*DY - XU2Y*DX;
    
    gu1xU1N = DY*gu1xXU1X;
    gu1yU1N = -DX*gu1yXU1Y;
    gt1U1N = XU1X*gt1DY + DY*gt1XU1X - XU1Y*gt1DX - DX*gt1XU1Y;
    gt2U1N = XU1X*gt2DY - XU1Y*gt2DX;
    gu2xU2N = DY*gu2xXU2X;
    gu2yU2N = -DX*gu2yXU2Y;
    gt1U2N = XU2X*gt1DY - XU2Y*gt1DX;
    gt2U2N = XU2X*gt2DY + DY*gt2XU2X - XU2Y*gt2DX - DX*gt2XU2Y;
    
    % G(ResF1) = G(ResF1) + CI3*U1N + CI6*U2N; % integro UN*N1
    % G(ResF2) = G(ResF2) + CI6*U1N + CI3*U2N; % integro UN*N2
    
    gG(ResF1,VarU1X) = gG(ResF1,VarU1X) + CI3*gu1xU1N;
    gG(ResF1,VarU1Y) = gG(ResF1,VarU1Y) + CI3*gu1yU1N;
    gG(ResF1,VarU2X) = gG(ResF1,VarU2X) + CI6*gu2xU2N;
    gG(ResF1,VarU2Y) = gG(ResF1,VarU2Y) + CI6*gu2yU2N;
    gG(ResF2,VarU1X) = gG(ResF2,VarU1X) + CI6*gu1xU1N;
    gG(ResF2,VarU1Y) = gG(ResF2,VarU1Y) + CI6*gu1yU1N;
    gG(ResF2,VarU2X) = gG(ResF2,VarU2X) + CI3*gu2xU2N;
    gG(ResF2,VarU2Y) = gG(ResF2,VarU2Y) + CI3*gu2yU2N;
    
    gG(ResF1,VarT1) = gG(ResF1,VarT1) + CI3*gt1U1N + CI6*gt1U2N;
    gG(ResF1,VarT2) = gG(ResF1,VarT2) + CI3*gt2U1N + CI6*gt2U2N;
    gG(ResF2,VarT1) = gG(ResF2,VarT1) + CI6*gt1U1N + CI3*gt1U2N;
    gG(ResF2,VarT2) = gG(ResF2,VarT2) + CI6*gt2U1N + CI3*gt2U2N;
    
end

end
