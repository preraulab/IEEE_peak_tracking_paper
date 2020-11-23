%CLIMSCALE Rescale the color limits of an image to remove outliers with percentiles
%
%   clim=climscale(h, ptiles);
%
%   h: image handle (optional, otherwise h=gca)
%   ptiles: percentiles (optional, default [5 98])
function clim=climscale(h, ptiles)
if nargin==0
    h=gca;
end

if nargin<2
    ptiles=[5 98];
end

children=get(h, 'children');

for i=1:length(children)
    if isprop(children(i),'cdata')
        data=get(children(i), 'cdata');
        data=data(:);
        bad_inds=isnan(data) | isinf(data);
        
        try
            clim=prctile(data(~bad_inds), ptiles);
        catch
            %Handle giant images
            
           data2=data(randi(length(data), 1, min(100000,length(data))));
           clear data;
           bad_inds=isnan(data2) | isinf(data2);
           
           clim=prctile(data2(~bad_inds), ptiles);
        end
        clim(2)=clim(2)+1e-10;
        set(h,'clim',clim);
    end
end