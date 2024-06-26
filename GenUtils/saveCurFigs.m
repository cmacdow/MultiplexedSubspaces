function saveCurFigs(handles,format,fn_base,save_dir,fancyflag)

%Camden MacDowell 2018
% @Synposis Saves figures in handles as type and to a specified dir. If dir is empty,
% saves to the current directory. To combine into a single pdf, save figs as
% pdfs and then use function 'CombinedPDFSInFolder'
% 
% 
% @Inputs
% @handles (required) 
% list of figure handes where handles = [figurehandl1 figurehandle2, ...]
%
% @format (required)
% Save format (e.g. '-dpdf' (or '-pdf' if fancyflag ==1) or any other format specificed in print fnct
% Updated 8-2021 to accept cell array for format
%
% @fn_base (required) 
% Base for the file name including the file extension. 
% This way all generated figures have the same
% suffix (easier to combined). Prefix for each file is 1 through numel(handels). 
% 
% @Savedir (optional)
% Directory to save the figures. if empty (def) use the current directory. 
%
% @fancyflag (optional)
% If you have complex images that are having issues saving correctly using typically saveas
% or if you want to save high quality images for publication purposes, set fancyflag to 1
% this will use the 'export_fig' function which is slow but better for custom rendering


if nargin <4
    save_dir = pwd; 
end
if nargin<5
	fancyflag=0;
end

if ~iscell(format)
    format = {format};
end

	

startDir = pwd; %so you return user to start dir after saving figs
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
cd(save_dir);

for cur_f = 1:numel(handles)
    %Save off the figure
    fig = handles(cur_f);
    set(fig,'Units', 'centimeters');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])   
    filename = sprintf('%d_%s',cur_f,fn_base);
    for cur_format = 1:numel(format)
        if fancyflag ==1
            export_fig(filename,format{cur_format},fig,'-nocrop','-transparent','-painters')
        else
            %Save the figure, with same suffix and iterate numbers as the prefix
           print(fig,sprintf('%d_%s',cur_f,fn_base),format{cur_format},'-painters','-r600')       
        end
    end %format loop
end %figure loop

cd(startDir) %return user to start dir

end %function
 


