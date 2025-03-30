function[] =  determine_disorder_scores(num_sequences, cohesive_only, prefix);
%%% HWdV 2020/05/13
% Processes disorder profiles, formatted as:
%                           column1=res#, column2 = disorder score [0,1]
% Inputs: 
% number of sequences (int), assumes that sequences are numbered
% cohesive_only (lower case bool), whether to only consider one domain (for NupY)
% prefix (str), repeating filename pattern.
%
% Outputs:
% log_disorder[full/cohesive_only].txt : text file containing what sequence has 100% of its
% residues disordered, together with its average disorder score
% graph1: plot of disorder score vs. sequence for all sequences
% graph2: plot of %disordered residues vs. design # incl. avg % of
%           disordered residues
% ranked_disorder[full/cohesive_only].txt: text file containing a sorted 
%                                          list with sequence numbers and 
%                                          their corresponding avg. disorder
%                                          score       
% 




    close all

    disorder_percentages = zeros(num_sequences,1);
    
    if cohesive_only == false
        fid = fopen('log_disorder_fullseq.txt','w');
        figure(1)
        hold on
        title('Disorder profiles, full sequence')
        xlabel('AA')
        ylabel('Disorder score (norm.)')
        box on
        pbaspect([1 1 1])
        
    elseif cohesive_only == true
        fid = fopen('log_disorder_cohesive_only.txt','w');
        figure(1)
        hold on
        title('Disorder profiles, cohesive domain')
        xlabel('AA')
        ylabel('Disorder score (norm.)')
        box on
        pbaspect([1 1 1])
    end



    for i=1:num_sequences
        open_str = [ prefix num2str(i) '.spotds.cleaned'];
        file_open = load(open_str,'-ascii');
        num_ordered = 0 ; 
        num_disordered = 0 ;
        
        if cohesive_only == false 
            for j=1:length(file_open)
                if file_open(j,2) >= 0.5
                    num_disordered = num_disordered +1 ;
                else
                    num_ordered = num_ordered +1 ;
                end
            end
            plot(file_open(:,1),file_open(:,2),'Color',[0.65 0 0 0.1],'LineWidth',1.5) ;
            disorder_percentages(i) = (num_disordered)/(num_disordered+num_ordered)*100;
            if disorder_percentages(i) == 100
                fprintf(fid, '%f\t%f\n', i, mean(file_open(:,2)) ) ;
            end
            
        elseif cohesive_only == true
            for j=1:613
                if file_open(j,2) >= 0.5
                    num_disordered = num_disordered +1 ;
                else
                    num_ordered = num_ordered +1 ;
                end
            end    
            plot(file_open(1:613,1),file_open(1:613,2),'Color',[0.65 0 0 0.1],'LineWidth',1.5) ;
            disorder_percentages(i) = (num_disordered)/(num_disordered+num_ordered)*100;
            if disorder_percentages(i) == 100
                fprintf(fid, '%f\t%f\n', i, mean(file_open(1:613,2)) ) ;
            end
        end

        
        
    end
    fclose(fid);
    
    
    figure(2)
    h2 = plot(linspace(1,num_sequences,num_sequences),disorder_percentages) ;
    avg_perc = mean(disorder_percentages);
    xlabel('Design #') ;
    ylabel(' Disordered residues (%)') ;
    h2.LineWidth = 2 ;
    h2.Color = [0.65 0 0] ;
    ylim([70,100])
    titlestring1=[ prefix '1:' num2str(num_sequences)] ;
    titlestring2=['avg score:' num2str(avg_perc) ] ;
    title({titlestring1,titlestring2},'Interpreter','none');
    pbaspect([1 1 1]);
    
    if cohesive_only == false
            diso_scores = load('log_disorder_fullseq.txt','-ascii');
            [diso_scores(:,2),idx] = sort(diso_scores(:,2));
            diso_scores(:,1) = diso_scores(idx,1);
            fid = fopen('ranked_disorder_fullseq.txt','w');
            for i = 1:length(diso_scores)
                fprintf(fid,'%f\t%f\n',diso_scores(i,1),diso_scores(i,2));
            end
            fclose(fid);

     elseif cohesive_only == true
            diso_scores = load('log_disorder_cohesive_only.txt','-ascii');
            [diso_scores(:,2),idx] = sort(diso_scores(:,2));
            diso_scores(:,1) = diso_scores(idx,1);
            fid = fopen('ranked_disorder_cohesive_only.txt','w');
            for i = 1:length(diso_scores)
                fprintf(fid,'%f\t%f\n',diso_scores(i,1),diso_scores(i,2));
            end
            fclose(fid)   ;         
     end
end

