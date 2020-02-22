%% REINITIALISATION DES VARIABLES ET FIGURES
clc
clear all
close all

%% exemples
 %% Plaque (0,90,90,0)
    
        Plaque = [] ;
        Plaque.e = [10 10 10 10]*1e-3 
        Plaque.NCouches = 4
        Plaque.THETA = [0 90 90 0]*pi/180 ;
        Plaque.NElmts =[10 10 10 10]
        Plaque.EL = [1.72e+11 1.72e+11 1.72e+11 1.72e+11 ]
        Plaque.ET = [6.89e+09 6.89e+09 6.89e+09 6.89e+09]
        Plaque.EN = [6.89e+09 6.89e+09 6.89e+09  6.89e+09]
        Plaque.GLT = [3.45e+09 3.45e+09 3.45e+09 3.45e+09]
        Plaque.GTN = [2.75e+09 2.75e+09 2.75e+09 2.75e+09 ]
        Plaque.GNL = [3.45e+09 3.45e+09 3.45e+09 3.45e+09 ]
        Plaque.nuLT = [0.25 0.25 0.25 0.25]
        Plaque.nuTN = [0.25 0.25 0.25 0.25 ]
        Plaque.nuLN = [0.25 0.25 0.25 0.25 ]
        Plaque.rho = [2260 2260 2260 2260]
        Plaque.nuTL = Plaque.nuLT.*Plaque.ET./Plaque.EL
        Plaque.nuNL = Plaque.nuLN.*Plaque.EN./Plaque.EL
        Plaque.nuNT = Plaque.nuTN.*Plaque.EN./Plaque.ET
        fmin = 50;
        fmax =70000;

 %% Plaque (0,-45,90,45,45,90,-45,0)
    
        Plaque = [] ;
        Plaque.e = [10 10 10 10 10 10 10 10]*1e-3 
        Plaque.NCouches = 8
        Plaque.THETA =[0 -45 90 45 45 90 -45 0]*pi/180;
        Plaque.NElmts = [10 10 10 10 10 10 10 10];
        Plaque.EL = [1.72e+11 1.72e+11 1.72e+11 1.72e+11 1.72e+11 1.72e+11 1.72e+11 1.72e+11]
        Plaque.ET = [6.89e+09 6.89e+09 6.89e+09 6.89e+09 6.89e+09 6.89e+09 6.89e+09 6.89e+09]
        Plaque.EN = [6.89e+09 6.89e+09 6.89e+09  6.89e+09 6.89e+09 6.89e+09 6.89e+09 6.89e+09]
        Plaque.GLT = [3.45e+09 3.45e+09 3.45e+09 3.45e+09 3.45e+09 3.45e+09 3.45e+09 3.45e+09]
        Plaque.GTN = [2.75e+09 2.75e+09 2.75e+09 2.75e+09 2.75e+09 2.75e+09 2.75e+09 2.75e+09]
        Plaque.GNL = [3.45e+09 3.45e+09 3.45e+09 3.45e+09 3.45e+09 3.45e+09 3.45e+09 3.45e+09]
        Plaque.nuLT = [0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25]
        Plaque.nuTN = [0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25]
        Plaque.nuLN = [0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25]
        Plaque.rho = [2260 2260 2260 2260 2260 2260 2260 2260]
        Plaque.nuTL = Plaque.nuLT.*Plaque.ET./Plaque.EL
        Plaque.nuNL = Plaque.nuLN.*Plaque.EN./Plaque.EL
        Plaque.nuNT = Plaque.nuTN.*Plaque.EN./Plaque.ET
        fmin = 50;
        fmax=70000;
%% Calculs EF
%% CALCUL DES NOMBRES D'ONDE
% Parametres
nF = 100 ;
scaleF = 'lin' ; 
phi0 =0;% % angle de propagation 
nPhi =1; %  Nombre de points sur la diagramme polaire
gammaMax = .1 ; % Critère sur l'amortissement pour le tri
        
    % Premiers calculs
        switch scaleF
            case 'lin'
                F = linspace(fmin,fmax,nF) ;
            case 'log'
                F = logspace(log10(fmin),log10(fmax),nF) ;
        end
        PHI = phi0+(0:nPhi-1)*2*pi/nPhi ;
        
    % CALCUL DES NOMBRES D'ONDE
        out = Shorter(Plaque,F,PHI) ; % Calcul nb ondes
        out = SortBranches3(out) ; % tri des branches
        out = TriPropagativ3(out,gammaMax) ; % tri des nb donde propagtifs
    
    % Résultats
        Data.K = out.k ;  %tenseur dordre trois k(f,phi,mode) = K[f,phi,mode]
        U = out.U ;
    % Calculs
        h = sum(Plaque.e) ;
        Data.Kstar = h*out.k;
        Data.C = 2*pi*repmat(reshape(out.F,[nF 1 1]),[1 1 size(out.k,3)])./out.k ;
        C2 = (Data.C).^2 ;
     
%% Ajout valeurs Reissner
        iPhi = 1;
        ondes = 1;
        nF=100;
        h = sum(Plaque.e) ;
        Plaque.Kappa = pi/sqrt(12) ;
        gammaMax = .1 ; % Critère sur l'amortissemen pour le tri
        rei = Reissner(Plaque,F,PHI) ; % Calcul nb ondes
        rei = SortBranches3(rei) ; % tri des branches
        rei = TriPropagativ3(rei,gammaMax) ; % tri des nb donde propagtifs
        starK_Reissner = h*rei.k;
        C_Reissner = 2*pi*repmat(reshape(rei.F,[nF 1 1]),[1 1 size(rei.k,3)])./rei.k ;       
        creiss=real(squeeze(C_Reissner(:,iPhi,ondes)));      
        % Ajout valeurs Kirchhoff
     C_Bend_LK =real(sqrt(2*pi*F).*(rei.Matrices.D(1,1)/rei.Matrices.M)^.25) ;
         
%% Facteur de normalisation
res  = calculs_elasticite(Plaque,PHI);
Nodes=res.Nodes; 
rhobar=res.rhobar %moyenne de rho
Gbar=res.Gbar %moyenne de G_LN
cs=sqrt(Gbar/rhobar) %facteur de normalisation

%% Résultats numériques
ratiokl=F./C_Bend_LK*h;
ckl=transpose(C_Bend_LK);
ratio=Data.Kstar(:,1,1)/(2*pi);
fosdt=creiss/cs; %R-M
fem=Data.C(:,1,1)/cs; % E-F
kl=C_Bend_LK/cs; %K-L

%% Plaque (0,90,90,0) %voir calculs calculs_BG
c_bg=(6.25e-34*(8.9453e134*Data.K(:,1,1).^12 + 3.2287e140*Data.K(:,1,1).^10 + 3.0881e145*Data.K(:,1,1).^8 + 3.1846e149*Data.K(:,1,1).^6 + 8.8962e152*Data.K(:,1,1).^4 + 1.0538e155*Data.K(:,1,1).^2 - 1.3902e108).^(1/2))./(1.7291e31*Data.K(:,1,1).^6 + 3.2008e36*Data.K(:,1,1).^4 + 1.7706e40*Data.K(:,1,1).^2 + 2.1443e42);

c_scp=(707.1067812*Data.K(:,1,1)).*sqrt((3.416964241*10^14+3.173975576*10^12*Data.K(:,1,1).^2+6.73070307*10^8*Data.K(:,1,1).^4).*(6.117953519*10^12+1.328420401*10^9*Data.K(:,1,1).^2))./(3.416964241*10^14+3.173975576*10^12*Data.K(:,1,1).^2+6.73070307*10^8*Data.K(:,1,1).^4);
c_ssp=(1414.213562*Data.K(:,1,1)).*sqrt((5.353200980*10^14+4.188474679*10^12*Data.K(:,1,1).^2+6.73070307*10^8*Data.K(:,1,1).^4).*(2.396179215*10^12+3.93306847*10^8*Data.K(:,1,1).^2))./(5.353200980*10^14+4.188474679*10^12*Data.K(:,1,1).^2+6.73070307*10^8*Data.K(:,1,1).^4);

bg=c_bg/cs;
scp=c_scp/cs;
ssp=c_ssp/cs;

%% Plaque (0,-45,90,45,45,90,-45,0)  %voir calculs calculs_BG
c_bg=  (0.5*(5.763e70*Data.K(:,1,1).^12 + 1.5003e75*Data.K(:,1,1).^10 + 1.018e79*Data.K(:,1,1).^8 + 5.421e81*Data.K(:,1,1).^6 + 8.5215e83*Data.K(:,1,1).^4 + 2.7304e85*Data.K(:,1,1).^2 + 4.1326e77).^(1/2))./(1.0822e32*Data.K(:,1,1).^6 + 1.4169e36*Data.K(:,1,1).^4 + 4.3021e38*Data.K(:,1,1).^2 + 1.5805e40);


c_scp=(1069.044968*Data.K(:,1,1)).*sqrt((8.174229540*10^12+2.335314781*10^11*Data.K(:,1,1).^2+8.13488071*10^8*Data.K(:,1,1).^4).*(1.954264779*10^11+8.12194499*10^8*Data.K(:,1,1).^2))./(8.174229540*10^12+2.335314781*10^11*Data.K(:,1,1).^2+8.13488071*10^8*Data.K(:,1,1).^4);
c_ssp=(-10000*Data.K(:,1,1)).*sqrt(-(1.*(-9.379487220*10^13-2.315477441*10^12*Data.K(:,1,1).^2-5.694416497*10^9*Data.K(:,1,1).^4)).*(2.562758123*10^10+7.1595949*10^7*Data.K(:,1,1).^2))./(-9.379487220*10^13-2.315477441*10^12*Data.K(:,1,1).^2-5.694416497*10^9*Data.K(:,1,1).^4);

bg=c_bg/cs;
scp=c_scp/cs;
ssp=c_ssp/cs;
%% Representation Graphique Des Resultats
set(gca,'FontSize',25)
pkl=plot(ratiokl,kl,'g','LineWidth',2); %resultats de KL
hold on
pfem=plot(ratio,fem,'k','LineWidth',2); %EF
hold on 
pfosdt=plot(ratio,fosdt,'b','LineWidth',2);%RM
hold on 
pbg=plot(ratio,bg,'r','LineWidth',2);%BG
hold on
pssp=plot(ratio,ssp,'--','Color',[.77 .02 .28],'LineWidth',2); %SSP
hold on
pscp=plot(ratio,scp,'--','Color',[.22 .45 .2],'LineWidth',2); %SCP
hold off
ylim([0 1])
xlim([0 0.5])
xlabel(' Rapport \''{e}paisseur-longueur d''onde $(h/\lambda)$','Interpreter','latex','FontSize',20);
ylabel('$\frac{\mbox{Vitesse des ondes de flexion}}{\mbox{Facteur de normalisation}}=\frac{c}{c_S}$','Interpreter','latex','FontSize',20);
lgd={'KL','FEM','$\frac{\pi^2}{12}$ FOSDT','BG','SSP','SCP'};
legendflex([pkl,pfem,pfosdt,pbg,pssp,pscp],lgd,'ref', gcf ,'anchor', {'s','s'},'buffer',[+80 +80], 'ncol',6,'fontsize',25);