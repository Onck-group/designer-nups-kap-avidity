%calc_spacing.m







[glfgpart,collapsed_no_GLFG]=regexp(nup116domain,'GLFG','match','split');



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

mean(spacings)
std(spacings)

sum_spacing/counter


%collapsed_no_GLFG = strjoin(collapsed_no_GLFG,'');

%[fgpart,collapsed_no_FGs]=regexp(collapsed_no_GLFG,'FG','match','split');

%collapsed_no_FGs=strjoin(collapsed_no_FGs,'');