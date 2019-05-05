function erps = plot_multi_ch( epochs, t, ys )

hold off;
if nargin < 3
    ys = 250;
end; if nargin == 1
    t = 1:size(epochs,2);
end

% hold off;
chanlist = 1:size(epochs, 1);
y_center = linspace( -ys, ys, length(chanlist) );
c=colormap('jet');
c_space = c( round(linspace(1,size(c,1),length(chanlist))), :) ;

% subplot(1,2,1);
erps = [];
for chanIdx = chanlist
    d = epochs( chanIdx, :) ...
        - (y_center(find(chanlist==chanIdx)) + nanmean( epochs( chanIdx, :)));
    plot( t, d, 'Color', c_space( find(chanlist==chanIdx),: )...
        , 'LineWidth', 1);
    if chanIdx == chanlist(1), hold on; end
end
set(gca, 'YTick', [] , 'Clipping', 'on');



function [idx,short] = hb_findIdx( range, fullData )
% function idx = hb_findIdx( range, fullData )

short=false;
idx = max(find( fullData < range(1)))+1:max(find(fullData<range(2)));

if length(idx)==0
%     idx
    idx = 1:max(find(fullData<range(2)));
end

if range(2) > fullData(end)
    short=1;
%     disp(['Requested range: .. ' num2str( range(2)) ]);
%     disp(['Actual range: .. ' num2str(fullData(end)) ]);
%     disp(['Short: ' num2str( short) ]); 
end

return





function str= hb_num2str( num )

if num > 9
    str = num2str( num );
else
    str = ['0' num2str( num )];
end
return

