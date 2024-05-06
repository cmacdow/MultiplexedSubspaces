function [neu_range,regD,locD,ratD,errinfo] = load_local_interregional_dim_data(filedir,method,thresh,normtype)
if nargin <2; method = 1; end %how to average the dimesnions across subsamples. def 1 = mean
if nargin <3; thresh = 0.8; end
if nargin <4; normtype = 'mean'; end

if ispc
    folder = ['Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\Analysis\',filedir,'\'];
else
    folder = ['/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/Analysis/',filedir,'\'];
end

regD = NaN(6,14,8,100); %interregional reliable dimensions
locD= NaN(6,14,8,100);  %local reliable dimensions
ratD = NaN(6,14,8,100); %ratio

errinfo=[]; %store runs that didn't finish

%load it into a nice lil' package 
for cur_rec = 1:6
    for cur_motif = 1:14
        for cur_area = 1:8
            if cur_rec == 3 && cur_area == 4 || cur_rec == 3 && cur_area == 6 || cur_rec ==4 && cur_area ==4
                 %missing areas on some recordings
            else
                try %to catch any missing recs
                fn = [folder,sprintf('rec%d_motif%d_area%d_thresh%g_norm%s.mat',cur_rec,cur_motif,cur_area,thresh,normtype)];                
                x = load(fn,'x_dim','x_dim_local','curve_dim','curve_dim_local','cur_rec','cur_motif','cur_area','n_subsample','num_neu','rPEV','rPEV_local','neu_range');
                
                %get the average across subsamples 
                d = x.curve_dim;
                dLoc = x.curve_dim_local;

                %get number of reliable dimensions
                d = getLastD(d);
                dLoc = getLastD(dLoc);

                %remove subsamples with no reliable dimensions
                d(d==0)=NaN;
                dLoc(dLoc==0)=NaN;

                if method ==1 %mean
                    regD(cur_rec,cur_motif,cur_area,1:size(x.x_dim,2)) = nanmean(d,2);
                    locD(cur_rec,cur_motif,cur_area,1:size(x.x_dim,2)) = nanmean(dLoc,2);
                    ratD(cur_rec,cur_motif,cur_area,1:size(x.x_dim,2)) = nanmean(dLoc./d,2);
                    
                elseif method ==0 %median
                    regD(cur_rec,cur_motif,cur_area,1:size(x.x_dim,2)) = nanmedian(d,2);
                    locD(cur_rec,cur_motif,cur_area,1:size(x.x_dim,2)) = nanmedian(dLoc,2);
                    ratD(cur_rec,cur_motif,cur_area,1:size(x.x_dim,2)) = nanmedian(dLoc./d,2); 
                    
                else
                    error('Unknown method');
                end
                catch 
                    fprintf('\ncur_rec = %d; cur_motif = %d; cur_area = %d;\n',cur_rec,cur_motif,cur_area)
                    errinfo = cat(1,errinfo,[cur_rec,cur_motif,cur_area]);
                end
            end
        end
    end
end

neu_range = x.neu_range;



end %function end




