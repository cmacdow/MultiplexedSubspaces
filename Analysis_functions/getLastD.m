function [Y,Yauc,X] = getLastD(X)
[x,~,z] = size(X);
X(X<=0)=NaN;

Y = NaN([x,z]);
Yauc = NaN([x,z]);
for xx = 1:x
    for zz = 1:z        
        if sum(~isnan(squeeze(X(xx,:,zz))))>0
            y = squeeze(X(xx,:,zz));
            y(find(isnan(y),1,'first'):end)=[];                        
            Y(xx,zz) = numel(y);
            %get AUC         
            if numel(y)>1
                y = [0,cumsum(y)/sum(y)];
                Yauc(xx,zz)= trapz(y/numel(y));
            end
        end
    end
end

end %subfunction end


%     y = [0,y(~isnan(y)),1];
%     %edge case, were - because of the cross validation the model performs
%     %better than the full model
%     y(y>1)=1;
%     auc = trapz(y/numel(y));


%old version
% for xx = 1:x
%     for yy = 1:y
%         for zz = 1:z
%             try %while troubleshooting
%             Y(xx,yy,zz) = find(~isnan(squeeze(X(xx,yy,zz,:))),1,"last");
%             catch
%             end
%         end
%     end
% end