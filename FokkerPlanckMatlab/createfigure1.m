function createfigure1(X1, YMatrix1)
%CREATEFIGURE1(X1,YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 16-Aug-2012 01:31:03

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'FontSize',15);
% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 7000000]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0 6e-007]);
box(axes1,'on');
hold(axes1,'all');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'Parent',axes1);
set(plot1(1),'Marker','.','DisplayName','maxwellian');
set(plot1(2),'Color',[0 0 1],'DisplayName','ne=1E16, n0=1E17');
set(plot1(3),'DisplayName','ne=1E16, n0=5E17');
set(plot1(4),'DisplayName','ne=1E16, n0=1E18','Color',[0 0 0]);
set(plot1(5),'LineStyle','--','Color',[0 0 1],...
    'DisplayName','ne=5E16, n0=1E17');
set(plot1(6),'LineStyle','--','Color',[1 0 0],...
    'DisplayName','ne=5E16, n0=5E17');
set(plot1(7),'LineStyle','--','DisplayName','ne=5E16, n0=1E18',...
    'Color',[0 0 0]);
set(plot1(8),'LineStyle',':','DisplayName','ne=1E17, n0=1E17');
set(plot1(9),'LineStyle',':','Color',[1 0 0],...
    'DisplayName','ne=1E17, n0=5E17');
set(plot1(10),'LineStyle',':','DisplayName','ne=1E17, n0=1E18',...
    'Color',[0 0 0]);

% Create xlabel
xlabel('v [m/s]','FontWeight','bold','FontSize',15);

% Create ylabel
ylabel('f [1/v]','FontWeight','bold','FontSize',15);

% Create title
title('Electron Velocity Distribution at T_{eff} = 7.5 eV',...
    'FontWeight','bold',...
    'FontSize',18);

% Create legend
legend(axes1,'show');
