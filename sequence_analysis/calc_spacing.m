%calc_spacing.m
clear all
close all


load('workspace_sequences.mat');



%lst_dir = dir('native-yeast-sequences/one-line-sequences/*domain.txt');

FG_Nup_spacings=zeros(length(FG_Nup_CHvals),1);
FG_Nup_spacings_std=zeros(length(FG_Nup_CHvals),1);

for nupnum=1:length(FG_Nup_names)
    filepathname=['native-yeast-sequences/one-line-sequences/' char(FG_Nup_names{nupnum}) '-domain.txt'];
    [mean_spacing,std_spacing] = calc_FG_spacing(filepathname);
    FG_Nup_spacings(nupnum)=mean_spacing;
    FG_Nup_spacings_std(nupnum)=std_spacing;
    
end


close all
figure()

hold all
box on
xlabel('C/H ratio (a.u.)');
ylabel('FG-spacing (residues)');

s1=scatter(NupY_CH_vals_calc,NupY_dFG_vals,30,'filled');
s1.MarkerEdgeColor=[0 0 0];
s2=scatter(FG_Nup_CHvals,FG_Nup_spacings,30,'filled');
s2.MarkerEdgeColor=[0 0 0];
text(FG_Nup_CHvals,FG_Nup_spacings,FG_Nup_names);
ylim([0,120]);

ll=line([0.0478 0.0478], [0 200]);
ll.LineStyle='--';
ll.Color=[0 0 0 0.8];

ll2=line([0 1], [13 13]);
ll2.LineStyle='--';
ll2.Color=[0 0 0 0.8];

figure()

hold all
box on
xlabel('C/H ratio (a.u.)');
ylabel('FG-spacing (residues)');

s1=scatter(NupY_CH_vals_calc,NupY_dFG_vals,30,'filled');
s1.MarkerEdgeColor=[0 0 0];
s2=scatter(FG_Nup_CHvals,FG_Nup_spacings,30,'filled');
s2.MarkerEdgeColor=[0 0 0];
text(FG_Nup_CHvals,FG_Nup_spacings,FG_Nup_names);

ll=line([0.0478 0.0478], [0 200]);
ll.LineStyle='--';
ll.Color=[0 0 0 0.8];

ll2=line([0 1], [13 13]);
ll2.LineStyle='--';
ll2.Color=[0 0 0 0.8];






function [mean_spacing,std_spacing] = calc_FG_spacing(filepathname) 

    
    %loadstr=[lst_dir(nup_index).folder '/' lst_dir(nup_index).name];
    loadstr=filepathname;
    fid=fopen(loadstr,'rt');
    nup_sequence=textscan(fid,'%s');
    %nup_sequence=load(loadstr,'-ascii');
    nup_sequence=char(nup_sequence{1});

    [glfgpart,collapsed_no_GLFG]=regexp(nup_sequence,'GLFG','match','split');



    [fgpart,collapsed_no_FG]=regexp(collapsed_no_GLFG,'FG','match','split');



    tot_cells=length(collapsed_no_FG);

    spacings = [] ;
    counter = 0 ;
    sum_spacing = 0;

    for i = 1:tot_cells

        len_cell = length(collapsed_no_FG{1,i});

        for j = 1:len_cell

            length_spacer = length(collapsed_no_FG{1,i}{j});

            spacings = [spacings ; length_spacer ] ;
            if length_spacer ~= 0 
                counter= counter+1 ;
                sum_spacing = sum_spacing + length_spacer ;

            end



        end

    end

    spacings(spacings == 0) = [] ;

    mean_spacing=mean(spacings);
    std_spacing=std(spacings);

    sum_spacing/counter;

end