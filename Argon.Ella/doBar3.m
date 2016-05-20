function doBar3(ymatrix1, level)
%CREATEFIGURE(YMATRIX1)
%  YMATRIX1:  bar matrix data

%  Auto-generated by MATLAB on 06-Oct-2011 16:12:41

% Create figure
figure1 = figure('Name','Contributions');
colormap('autumn');

% Create axes
axes1 = axes('Parent',figure1,...
    'XTickLabel',{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25'},...
    'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25],...
    'XGrid','on');
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 25]);
box(axes1,'on');
hold(axes1,'all');

% Create multiple lines using matrix input to bar
bar1 = bar(ymatrix1,'Parent',axes1);
set(bar1(1),'DisplayName','n_e = 10^6 cm^{-3}');
set(bar1(2),'DisplayName','n_e = 10^{11} cm^{-3}');

% Create xlabel
xlabel('Contributing Level');

% Create ylabel
ylabel('Relative Contribution');

% Create title
title(['Level ' num2str(level)]);

% Create legend
legend(axes1,'show');
