function [b,area,area_label,a] = loadVisBeta(data,cur_d,cur_rec,cur_m,targ_area)
%see PlotBetasFigure and CompareEphysCalciumNetowrks
[b, area_label] = LoadVariable(data(cur_rec),'rrr_beta',targ_area,cur_d); %load betas
area_label = area_label(strcmp(area_label,targ_area)==0);
area=LoadVariable(data(cur_rec),'beta_region',targ_area);
%get one rec and one motif
b = (squeeze(b(cur_m,:)));
area = squeeze(area(cur_m,:));
idx = isnan(area);
area(idx)=[];
b(idx)=[];

%reorder by decreasing median beta
a = b;
[b,reord] = sort(b,'descend');
area = area(reord);
end