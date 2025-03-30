function [] = save_plot(outname,filetype)



    set(gcf,'renderer','painters')
    %print -depsc2 -tiff -r600 velmap_5_kaps.eps
    if strcmp(filetype,'png')==1
        print(outname,'-dpng','-r600')        
    end
    
    if strcmp(filetype,'eps')==1
        print(outname,'-depsc2','-tiff','-r300')        
    end

end