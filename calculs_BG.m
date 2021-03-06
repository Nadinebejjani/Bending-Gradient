close all
clear all

%% SYMBOLES UTILISES
syms a b c d k W w

%% Plaque (0,90,90,0) (dim S=4)

N=[-1.372e-08       0.8093   -6.226e-08      0.58409   -0.0015497    -0.062187;
  -9.7044e-09    0.0013389   7.5534e-08     -0.10774    -0.024767     -0.99387;
      0.46772   1.7097e-08      0.88373   8.7148e-08     0.016086  -0.00040081;
   -0.0025312   4.3046e-11     0.019544   2.7442e-09      -0.9995     0.024908;
      0.88387   1.9681e-08     -0.46759  -6.0345e-08    -0.011374   0.00028341;
  -2.4303e-08      0.58739   7.5552e-08     -0.80451    0.0021916     0.087946]

H=[-3.9306e+07        71884   3.0349e+08  -5.7842e+06  -1.5345e+10   2.3597e+09;
   7.2501e+06        118.8  -5.5979e+07      -9577.2   2.8319e+09  -4.2542e+08;
   5.0223e+08        18.36   9.4895e+08       93.456   1.7084e+07   -2.567e+06;
   1.1107e+07      0.40601   2.0986e+07       2.0692   3.7779e+05       -56767;
  -2.6573e+08      -9.7118  -5.0209e+08      -49.653  -9.0388e+06   1.3582e+06;
   5.4138e+07        52172  -4.1801e+08  -4.1983e+06   2.1153e+10  -3.1391e+09]

Ds=[8.0929e+05       9209.7  -9.9269e-14;
       9209.7   1.4719e+05  -9.4565e-12;
  -9.9269e-14  -9.4565e-12        36800]

rhobar=90.4

 Phiet= [ a; b; c; d; 0; 0 ]

%% Plaque (0,-45,90,45,45,90,-45,0)

%matrice H   
H=[969293326909.75269        2517811685142.3931       -942753252212.51624        68552438459.070755        665005932115.27979       -1368384455311.8342 ;    
   2517811685142.3931        6868377468993.8604       -2281480252761.4038       -1765775058149.8875        1613320213069.3728       -3556801775772.9941  ;   
  -942753252212.51611       -2281480252761.4038        1008467124077.5840       -1112972095045.3022       -709048142205.90735        1330556480403.0298  ;   
   68552438459.070717       -1765775058149.8875       -1112972095045.3022        12413638460720.896        761483971026.90759       -88671206937.805984 ;    
   665005932115.27979        1613320213069.3728       -709048142205.90735        761483971026.90747        498740196139.47717       -938557218741.70093 ;    
  -1368384455311.8342       -3556801775772.9946        1330556480403.0298       -88671206937.805847       -938557218741.70093        1932243699220.0830  ]
    
Ds=[4940858.2173107853        613962.30387853389       -702260.56560566509   ;  
   613962.30387853389        1630370.0637659924       -702260.67248213512  ;   
  -702260.56560566509       -702260.67248213512        1374969.0166635716   ]

rhobar=180.8
 
 %% CALCULS BG
%% Plaque (0,90,90,0) (dim S=4)
%Phi
Phi=vpa(N*Phiet,5)

%igradW
igradW =vpa([-(1i)*k*W; 0 ; 0 ; 0 ;0; -(sqrt(2)/2)*(1i)*k*W],5)

%Gamma
Gamma=vpa(Phi+igradW,5)

% tenseur R
R=vpa(H*Gamma,5)

% tenseur chi
Chi=vpa([-i*k*Phi(1);-i*k*Phi(2);-i*k*Phi(3)],5)

% tenseur M et gradM
M=vpa(Ds*Chi,5)

gradM=vpa([-i*k*M(1);-i*k*M(2);-i*k*M(3);0;0;0],5)

% Equations
E=vpa(R-gradM,5)

%simplifications 
EQ1=E(1)+(1/sqrt(2))*E(6)
EQ2=E(2)
EQ3=E(5)+(1/sqrt(2))*E(3)
EQ4=E(4)


Q1d1=-i*k*R(1)-i*k*R(6)/sqrt(2)

EQ5= w^2*rhobar*W+Q1d1

%equations to matrix form
eqns = [EQ5 == 0
        EQ1 == 0
        EQ2 == 0
        EQ3 == 0
        EQ4 == 0];
vars = [W a b c d];
[Mat,Vect] = equationsToMatrix(eqns,vars)  
 
% MAtrice A et Determinant de A
A=vpa(Mat,5)
detA=vpa(det(A),5)

% Resoudre l'equation detA=0 % w en fct de k
eqn = detA ==0
soleqn = solve(eqn,w)

wsol=vpa(soleqn,5)
%% Plaque (0,-45,90,45,45,90,-45,0)
A11=w^2*rhobar-k^2*(H(1,1)+sqrt(2)*H(6,1)+0.5*H(6,6));
A12=-i*k*(H(1,1)+(sqrt(2)/2)*H(1,6));
A13=-i*k*(H(1,2)+(sqrt(2)/2)*H(2,6));
A14=-i*k*(H(1,3)+(sqrt(2)/2)*H(3,6));
A15=-i*k*(H(1,4)+(sqrt(2)/2)*H(4,6));
A16=-i*k*(H(1,5)+(sqrt(2)/2)*H(5,6));
A17=-i*k*(H(1,6)+(sqrt(2)/2)*H(6,6));

A22=H(1,1)+k^2*Ds(1,1);
A23=H(1,2)+k^2*Ds(1,2);
A24=H(1,3)+k^2*Ds(1,3);
A25=H(1,4);
A26=H(1,5);
A27=H(1,6);

A33=H(2,2)+k^2*Ds(2,2);
A34=H(2,3)+k^2*Ds(2,3);
A35=H(2,4);
A36=H(2,5);
A37=H(2,6);

A44=H(3,3)+k^2*Ds(3,3);
A45=H(3,4);
A46=H(3,5);
A47=H(3,6);

A55=H(4,4);
A56=H(4,5);
A57=H(4,6);

A66=H(5,5);
A67=H(5,6);

A77=H(6,6); 

Mat=[ A11 A12 A13 A14 A15 A16 A17;
    A12 A22 A23 A24 A25 A26 A27;
    A13 A23 A33 A34 A35 A36 A37;
    A14 A24 A34 A44 A45 A46 A47;
    A15 A25 A35 A45 A55 A56 A57;
    A16 A26 A36 A46 A56 A66 A67;
    A17 A27 A37 A47 A57 A67 A77]

% Matrice A et Determinant de A
A=vpa(Mat,5)
detA=vpa(det(A),5)

% Resoudre l'equation detA=0
eqn = detA ==0
soleqn = solve(eqn,w)

wsol=vpa(soleqn,5)
