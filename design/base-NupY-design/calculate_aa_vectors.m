%calculate_aa_vectors.m


%SET AS HYDROPHOBIC: A,I,L,M,F,V,P,G or: AILFWV acc. to yamada (indeed, M,
%P, G are a bit odd)

%read in all of the different Nup domains.
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


%calculate the average c/h ratios, charged fractions ,hydrophobic
%fractions)
collapsedDomainHydro=[hphobicNup49C,hphobicNup57C,hphobicNup100C,hphobicNup116C,hphobicNup145NC];
collapsedDomainCharge=[chargedNup49C,chargedNup57C,chargedNup100C,chargedNup116C,chargedNup145NC];

extendedDomainHydro=[hphobicNup100E,hphobicNup116E,hphobicNup145NE];
extendedDomainCharge=[chargedNup100E,chargedNup116E,chargedNup145NE];

avgCDH=mean([hphobicNup49C,hphobicNup57C,hphobicNup100C,hphobicNup116C,hphobicNup145NC]);
avgCDC=mean([chargedNup49C,chargedNup57C,chargedNup100C,chargedNup116C,chargedNup145NC]);

avgEDH=mean([hphobicNup100E,hphobicNup116E,hphobicNup145NE]);
avgEDC=mean([chargedNup100E,chargedNup116E,chargedNup145NE]);

avgCHratioC=mean(avgCDC/avgCDH);
avgCHratioE=mean(avgEDC/avgEDH);


%initialize the AA vector
aaVector = ['R'   'K'     'D'    'E'    'N' 'Q' 'S' 'T' 'H' 'C' 'Y' 'G' 'A' 'P' 'M' 'V' 'W' 'I' 'L' 'F'];
%concatenate all the extended, collapsed domains

collapsed_domains_AA_pool  = strcat(Nup100C, Nup116C, Nup145NC, Nup49C, Nup57C) ;
extended_domains_AA_pool = strcat(Nup100E, Nup116E, Nup145NE)  ;


%correct for GLFG, FG
[glfgpart,collapsed_no_GLFG]=regexp(collapsed_domains_AA_pool,'GLFG','match','split');
collapsed_no_GLFG = strjoin(collapsed_no_GLFG,'');
[fgpart,collapsed_no_FGs]=regexp(collapsed_no_GLFG,'FG','match','split');

collapsed_no_FGs=strjoin(collapsed_no_FGs,'');



composition_AA_extended  = zeros(length(aaVector));
composition_AA_collapsed = zeros(length(aaVector));
composition_AA_collapsed_noFG = zeros(length(aaVector));

for i = 1:length(aaVector)
    currentAA=aaVector(i) ; 
    composition_AA_extended(i)  = count(extended_domains_AA_pool, currentAA)  ;
    composition_AA_collapsed(i) = count(collapsed_domains_AA_pool, currentAA) ;
    composition_AA_collapsed_noFG(i) = count(collapsed_no_FGs, currentAA) ;
end
composition_AA_extended  = composition_AA_extended/length(extended_domains_AA_pool)   ;
composition_AA_collapsed = composition_AA_collapsed/length(collapsed_domains_AA_pool) ;
composition_AA_collapsed_noFG = composition_AA_collapsed_noFG/length(collapsed_no_FGs) ;

fid = fopen('output_GLFG_Nup_sequence_analysis.dat','w');

fprintf(fid,'AA\tC_dom\tnoFG\tE_dom\n') ;
for i = 1:length(aaVector)
    fprintf(fid,'%s\t%f\t%f\t%f\n', aaVector(i), composition_AA_collapsed(i),composition_AA_collapsed_noFG(i),composition_AA_extended(i));
end


fclose(fid) ;



