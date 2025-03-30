%plot_scatter_properties.m

close all

            

fid=fopen('outcomes_kd_sorted.dat','r');
master_array=textscan(fid,'%f%f%f%f%f%f%f%f','Headerlines',0,'delimiter','\t');
fclose(fid);
%master_array=load('dG_values.dat','-ascii');
dFG_vals=master_array{1};
CH_vals=master_array{2};
kd_vals=master_array{3};    %use horzcat()

color_array={[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840]};




%type 1
figure()
us = scatter(CH_vals,dFG_vals,50,kd_vals,'Filled');
cc = colorbar();
cc.Label.String='{K_D} ({\mu}M)';
cc.Ruler.Scale = 'log';
cc.Ruler.MinorTick = 'on';
box on
xlabel('C/H ratio')
ylabel('FG-spacing')
pbaspect([1 1 1])





figure()

CH_iso = [0.024, 0.048, 0.096, 0.191, 0.382, 0.574] ;
ms=horzcat(CH_vals,dFG_vals,kd_vals);

for jj = 1:length(CH_iso)
    
    
    selector= ms(:,1) == CH_iso(jj) ;
    dfg=ms(selector,2);
    kd_p=ms(selector,3);
    subplot(1,2,1)
    hold all
    ss1=scatter(dfg,kd_p,50,'Filled');
    ss1.MarkerFaceColor=color_array{jj};
    ss1.MarkerEdgeColor=[0 0 0];
    pp=plot(dfg,kd_p);
    pp.Color=color_array{jj};
    
    ylabel('{K_D} ({\mu}M)');
    xlabel('FG-spacing (residues)')
    box on
    set(gca,'yscale','log');
    ylim([0,1000]);


    
end

%slegend('C/H=0.024', 'C/H=0.048', 'C/H=0.096', 'C/H=0.191', 'C/H=0.382', 'C/H=0.574', 'Location','SouthEast');





dFG_iso = [7, 13, 26, 52, 104];

for jj =1:length(dFG_iso)
    selector= ms(:,2) == dFG_iso(jj) ;
    dch=ms(selector,1);
    kd_p=ms(selector,3);
    subplot(1,2,2)

    hold all

    ss2=scatter(dch,kd_p,50,'Filled');
    ss2.MarkerFaceColor=color_array{jj};
    ss2.MarkerEdgeColor=[0 0 0];
    pp=plot(dch,kd_p);
    pp.Color=color_array{jj};
    ylim([0,1000]);
    set(gca,'yscale','log');
    ylabel('{K_D} ({\mu}M)');
    xlabel('C/H ratio (dimensionless)');
    box on
    
end

%legend('{d_{FG}}=7','{d_{FG}}=13','{d_{FG}}=26','{d_{FG}}=52','{d_{FG}}=104', 'Location','NorthEast');
