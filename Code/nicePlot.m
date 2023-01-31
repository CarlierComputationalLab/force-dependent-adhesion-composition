% Function that takes a matlab plot and makes it nice. If you want to save
% it, the argument 'save' should be 1. 

function nicePlot(figureHandle, fontSize, picturewidth, hw_ratio, saveDirectory, fname, save, showGraphs)
hfig = figureHandle;
hfig.Visible = showGraphs;
% Set the directory where the figure needs to be saved. fname is the name of the figure
fname = fullfile(saveDirectory,fname);
% picturewidth = 20; % set this parameter and keep it forever
% hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',fontSize) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
% print(hfig,fname,'-dpdf','-painters','-fillpage')
if save == 1
    print(hfig,fname, '-dsvg','-r350','-painters')
    print(hfig,fname, '-djpeg','-painters')
end

end
