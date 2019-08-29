%% DEMONSTRATION DU FONCTIONNEMENT DES SCRIPTS POUR LE CALCUL DE LA PROPAGATION DES ONDES


%% REINITIALISATION DES VARIABLES ET FIGURES

clc
clear all
close all

 %%
file = 'Plaque1.mat' ;
Data = load(file) ;

 %%
file = 'Plaque2.mat' ;
Data = load(file) ;
 %%
file = 'Plaque3.mat' ;
Data = load(file) ;
 %%
file = 'Plaque4.mat' ;
Data = load(file) ;

 %% Ajout valeurs Reissner
        iPhi = 1;
        ondes = 1;
        nF=1000;
       h = sum(Data.Plaque.e) ;
        Data.Plaque.Kappa = pi/sqrt(12) ;
        gammaMax = .1 ; % Critère sur l'amortissement pour le tri
        %Plaque.Kappa = sqrt(5/6) ;
        rei = Reissner(Data.Plaque,Data.F,Data.PHI) ; % Calcul nb ondes
        rei = SortBranches3(rei) ; % tri des branches
        rei = TriPropagativ3(rei,gammaMax) ; % tri des nb donde propagtifs
        %real(squeeze(C_Reissner(:,iPhi,ondes)))
        starK_Reissner = h*rei.k;
        C_Reissner = 2*pi*repmat(reshape(rei.F,[nF 1 1]),[1 1 size(rei.k,3)])./rei.k ;       
        creiss=real(squeeze(C_Reissner(:,iPhi,ondes)));      
        % Ajout valeurs Kirchhoff
     C_Bend_LK =real(sqrt(2*pi*Data.F).*(rei.Matrices.D(1,1)/rei.Matrices.M)^.25) ;
%% Calcul analytique

res  = analyticcalculus(Data.Plaque,Data.PHI);
intrhow=res.intrhow
intrho=res.intrho
Gmoy=res.Gmoy;
termkl=res.termkl
Gmoy=res.Gmoy
Gmoytest=res.Gmoytest 
%% figure des 5 modes ELEMENTS FINIS
%ondes de flexion
XData1=(Data.Kstar(:,1,1))/(2*pi);%ratio
YData1 = Data.C(:,1,1);%/cs;
%ondes de cisaillemen( in plane shear)
XData2=(Data.Kstar(:,1,2))/(2*pi);%ratio
YData2 = Data.C(:,1,2);%/cs;/cs;
%ondes de compression
XData3=(Data.Kstar(:,1,3))/(2*pi);%ratio
YData3 = Data.C(:,1,3);%/cs;/cs;
%ondes de cisaillement transverse
XData4=(Data.Kstar(:,1,4))/(2*pi);%ratio
YData4 = Data.C(:,1,4);%/cs;/cs;
%ondes 5
XData5=(Data.Kstar(:,1,5))/(2*pi);%ratio
YData5 = Data.C(:,1,5);%/cs;/cs;

pfem=plot(XData1,YData1,'k','LineWidth',2);
set(gca,'FontSize',25)
hold on
pord0=plot(XData2,YData2,'Color',[.5 .00 .5],'LineWidth',2);
hold on 
pord1=plot(XData3,YData3,'r','LineWidth',2);
hold on
pord2=plot(XData4,YData4,'b','LineWidth',2);
hold off
%xlim([0 0.5])
%ylim([0 10])
ylabel('$\frac{c}{c_S}$','Interpreter','latex','FontSize',20);
xlabel(' Thickness-wavelength ratio $(h/\lambda)$','Interpreter','latex','FontSize',20);%lgd={'FEM','p=0','p=1','p=2','p=3','p=4','p=5','p=6','p=7','p=8','p=9','p=10'};
lgd={'ondes 1','ondes 2','ondes 3','ondes 4'};
legendflex([pfem,pord0,pord1,pord2],lgd,'ref', gcf ,'anchor', {'n','n'},'buffer',[-80 -80], 'ncol',6,'fontsize',25);
set(gca,'xscale','log','yscale','log')
%% Definitions communes
cs=sqrt(Gmoy/intrhow); %Normalization factor
%KL
ckl=transpose(C_Bend_LK);
kl=ckl/cs;
ratiokl=Data.F./C_Bend_LK*h;
% axe des x
ratio=Data.Kstar(:,1,1)/(2*pi);
% EF
fem=Data.C(:,1,1)/cs;
%RM
fosdt=creiss/cs;


%% Plaque Une
cbg=(6.25e-34*(8.9453e134*Data.K(:,1,1).^12 + 3.2287e140*Data.K(:,1,1).^10 + 3.0881e145*Data.K(:,1,1).^8 + 3.1846e149*Data.K(:,1,1).^6 + 8.8962e152*Data.K(:,1,1).^4 + 1.0538e155*Data.K(:,1,1).^2 - 1.3902e108).^(1/2))./(1.7291e31*Data.K(:,1,1).^6 + 3.2008e36*Data.K(:,1,1).^4 + 1.7706e40*Data.K(:,1,1).^2 + 2.1443e42);
bg=cbg/cs;

%% Plaque Deux
cbg=  (2.5e-24*(6.4882e117*Data.K(:,1,1).^12 + 5.6277e123*Data.K(:,1,1).^10 + 1.2182e129*Data.K(:,1,1).^8 + 2.3835e132*Data.K(:,1,1).^6 + 1.3377e135*Data.K(:,1,1).^4 + 1.587e137*Data.K(:,1,1).^2 + 1.0807e99).^(1/2))./(1.8298e32*Data.K(:,1,1).^6 + 8.5673e37*Data.K(:,1,1).^4 + 9.4039e40*Data.K(:,1,1).^2 + 1.2924e43);
bg=cbg/cs;

%% Plaque Trois 
cbg=  (0.5*(5.763e70*Data.K(:,1,1).^12 + 1.5003e75*Data.K(:,1,1).^10 + 1.018e79*Data.K(:,1,1).^8 + 5.421e81*Data.K(:,1,1).^6 + 8.5215e83*Data.K(:,1,1).^4 + 2.7304e85*Data.K(:,1,1).^2 + 4.1326e77).^(1/2))./(1.0822e32*Data.K(:,1,1).^6 + 1.4169e36*Data.K(:,1,1).^4 + 4.3021e38*Data.K(:,1,1).^2 + 1.5805e40);
bg=cbg/cs;

%% Plaque Quatre
cbg=(2.5e-24*(1.1198e117*Data.K(:,1,1).^12 + 2.0545e121*Data.K(:,1,1).^10 + 1.0615e125*Data.K(:,1,1).^8 + 1.114e128*Data.K(:,1,1).^6 + 3.7531e130*Data.K(:,1,1).^4 + 2.7523e132*Data.K(:,1,1).^2 - 5.0049e94).^(1/2))./(6.6503e31*Data.K(:,1,1).^6 + 6.4123e35*Data.K(:,1,1).^4 + 4.0918e38*Data.K(:,1,1).^2 + 3.4713e40);
bg=cbg/cs;
%% Representation Graphique Des Resultats
set(gca,'FontSize',25)
%hold on
pkl=plot(ratiokl,kl,'g','LineWidth',2);
hold on
pfem=plot(ratio,fem,'k','LineWidth',2);
hold on 
pfosdt=plot(ratio,fosdt,'b','LineWidth',2);%'b-x'
hold on 
pbg=plot(ratio,bg,'r','LineWidth',2);
% hold on
% pssp=plot(ratio,ssp,'--','Color',[.77 .02 .28],'LineWidth',2);
% hold on
% pscp=plot(ratio,scp,'--','Color',[.22 .45 .2],'LineWidth',2);
hold off
ylim([0 1])
xlim([0 0.5])
 xlabel(' Thickness-wavelength ratio $(h/\lambda)$','Interpreter','latex','FontSize',20);
 ylabel('$\frac{\mbox{Flexural wave velocity}}{\mbox{Normalization factor}}=\frac{c}{c_S}$','Interpreter','latex','FontSize',20);
lgd={'KL','FEM','$\frac{\pi^2}{12}$ FOSDT','BG'};
legendflex([pkl,pfem,pfosdt,pbg],lgd,'ref', gcf ,'anchor', {'s','s'},'buffer',[+80 +80], 'ncol',6,'fontsize',25);

