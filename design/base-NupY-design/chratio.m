%chratiocalculator

%SET AS HYDROPHOBIC: A,I,L,M,F,V,P,G or: AILFWV acc. to yamada (indeed, M,
%P, G are a bit odd)

Nup49C=fileread('./Native_GLFG/Nup49coll');
hphobicNup49C=count(Nup49C,["A","I","L","F","W","V"])/(length(Nup49C));
chargedNup49C=count(Nup49C,["D","E","K","R"])/(length(Nup49C));
chRatioNup49C=chargedNup49C/hphobicNup49C;


Nup57C=fileread('./Native_GLFG/Nup57coll');
hphobicNup57C=count(Nup57C,["A","I","L","F","W","V"])/(length(Nup57C));
chargedNup57C=count(Nup57C,["D","E","K","R"])/(length(Nup57C));
chRatioNup57C=chargedNup57C/hphobicNup57C;



Nup116E=fileread('./Native_GLFG/Nup116ext');
Nup116C=fileread('./Native_GLFG/Nup116coll');
hphobicNup116E=count(Nup116E,["A","I","L","F","W","V"])/(length(Nup116E));
hphobicNup116C=count(Nup116C,["A","I","L","F","W","V"])/(length(Nup116C));
chargedNup116E=count(Nup116E,["D","E","K","R"])/(length(Nup116E));
chargedNup116C=count(Nup116C,["D","E","K","R"])/(length(Nup116C));
chRatioNup116E=chargedNup116E/hphobicNup116E;
chRatioNup116C=chargedNup116C/hphobicNup116C;


Nup100E=fileread('./Native_GLFG/Nup100ext');
Nup100C=fileread('./Native_GLFG/Nup100coll');
hphobicNup100E=count(Nup100E,["A","I","L","F","W","V"])/(length(Nup100E));
hphobicNup100C=count(Nup100C,["A","I","L","F","W","V"])/(length(Nup100C));
chargedNup100E=count(Nup100E,["D","E","K","R"])/(length(Nup100E));
chargedNup100C=count(Nup100C,["D","E","K","R"])/(length(Nup100C));
chRatioNup100E=chargedNup100E/hphobicNup100E;
chRatioNup100C=chargedNup100C/hphobicNup100C;



Nup145NE=fileread('./Native_GLFG/Nup145Next');
Nup145NC=fileread('./Native_GLFG/Nup145Ncoll');
hphobicNup145NE=count(Nup145NE,["A","I","L","F","W","V"])/(length(Nup145NE));
hphobicNup145NC=count(Nup145NC,["A","I","L","F","W","V"])/(length(Nup145NC));
chargedNup145NE=count(Nup145NE,["D","E","K","R"])/(length(Nup145NE));
chargedNup145NC=count(Nup145NC,["D","E","K","R"])/(length(Nup145NC));
chRatioNup145NE=chargedNup145NE/hphobicNup145NE;
chRatioNup145NC=chargedNup145NC/hphobicNup145NC;

NupXE=fileread('./Native_GLFG/NupXext');
NupXC=fileread('./Native_GLFG/NupXcoll');
hphobicNupXE=count(NupXE,["A","I","L","F","W","V"])/(length(NupXE));
hphobicNupXC=count(NupXC,["A","I","L","F","W","V"])/(length(NupXC));
chargedNupXE=count(NupXE,["D","E","K","R"])/(length(NupXE));
chargedNupXC=count(NupXC,["D","E","K","R"])/(length(NupXC));
chRatioNupXE=chargedNupXE/hphobicNupXE;
chRatioNupXC=chargedNupXC/hphobicNupXC;



collapsedDomainHydro=[hphobicNup49C,hphobicNup57C,hphobicNup100C,hphobicNup116C,hphobicNup145NC];
collapsedDomainCharge=[chargedNup49C,chargedNup57C,chargedNup100C,chargedNup116C,chargedNup145NC];

extendedDomainHydro=[hphobicNup100E,hphobicNup116E,hphobicNup145NE];
extendedDomainCharge=[chargedNup100E,chargedNup116E,chargedNup145NE];



NupYC=fileread('./NupYC_1.txt');
hphobicNupYC=count(NupYC,["A","I","L","F","W","V"])/(length(NupYC));
chargedNupYC=count(NupYC,["D","E","K","R"])/(length(NupYC));
chRatioNupYC=chargedNupYC/hphobicNupYC;

NupYE=fileread('./NupYE_1.txt');
hphobicNupYE=count(NupYE,["A","I","L","F","W","V"])/(length(NupYE));
chargedNupYE=count(NupYE,["D","E","K","R"])/(length(NupYE));
chRatioNupYE=chargedNupYE/hphobicNupYE;


avgCDH=[hphobicNup49C,hphobicNup57C,hphobicNup100C,hphobicNup116C,hphobicNup145NC];
avgCDC=[chargedNup49C,chargedNup57C,chargedNup100C,chargedNup116C,chargedNup145NC];



avgEDH=[hphobicNup100E,hphobicNup116E,hphobicNup145NE];
avgEDC=[chargedNup100E,chargedNup116E,chargedNup145NE];

avgCHratioC=mean(avgCDC/avgCDH);
avgCHratioE=mean(avgEDC/avgEDH);

fprintf('C/H ratio collapsed domains GLFG Nup avg:\t %f\n',avgCHratioC)
fprintf('C/H ratio collapsed domain NupY avg:\t %f\n',chRatioNupYC)
fprintf('relative error:\t%f%%\n',100*(chRatioNupYC-avgCHratioC)/avgCHratioC)

fprintf('C/H ratio extended domains GLFG Nup avg: %f\n',avgCHratioE)
fprintf('C/H ratio extended domain NupY avg: %f\n',chRatioNupYE)
fprintf('relative error:\t%f%%\n',100*(chRatioNupYE-avgCHratioE)/avgCHratioE)

hold on
%xlim([0,40])
pbaspect([1 1 1])
ylim([20,30])
set(gca,'box','on')
xlabel('Charged AA (%)') %DEKR
ylabel('Hydrophobic AA (%)') ; %AILFW


s1=scatter(extendedDomainCharge*100,extendedDomainHydro*100,'filled');


tl1=text(mean(extendedDomainCharge)*100,28.5,'Extended domains','HorizontalAlignment','center');

s2=scatter(collapsedDomainCharge*100,collapsedDomainHydro*100,'filled');

tl2=text(mean(collapsedDomainCharge)*100,25,'Collapsed domains','HorizontalAlignment','Left');


sNupYc=scatter(chargedNupYC*100,hphobicNupYC*100,'filled');

sNupYc.MarkerEdgeColor=[0 0 0];
sNupYc.CData=[0, 0, 0];

sNupYe=scatter(chargedNupYE*100,hphobicNupYE*100,'filled');

sNupYe.MarkerEdgeColor=[0 0 0];
sNupYe.CData=[0, 0, 0];


sNupXc=scatter(chargedNupYC*100,hphobicNupXC*100,'filled');

sNupXc.MarkerEdgeColor=[0 0 0];
sNupXc.CData=[0, 1, 0];

sNupXe=scatter(chargedNupYE*100,hphobicNupXE*100,'filled');

sNupXe.MarkerEdgeColor=[0 0 0];
sNupXe.CData=[1, 0, 0];



s1.CData=[1, 0.898, 0.906];
s1.MarkerEdgeColor=[0 0 0];

s2.CData=[0.8 1 0.816];
s2.MarkerEdgeColor=[0 0 0];

tl1.Color=[0.64, 0 0];
tl2.Color=[0 0.65 0];



