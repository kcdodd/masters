function plotRadialSolution(xdata1, ydata1, zdata1, U1, V1)
%CREATEFIGURE(XDATA1,YDATA1,ZDATA1,U1,V1)
%  XDATA1:  contour x
%  YDATA1:  contour y
%  ZDATA1:  contour z
%  U1:  quiver u
%  V1:  quiver v

%  Auto-generated by MATLAB on 29-May-2012 19:02:58

% Create figure
figure1 = figure('Colormap',...
    [0 0 0;0.00584795325994492 0.00584795325994492 0.00584795325994492;0.0116959065198898 0.0116959065198898 0.0116959065198898;0.0175438597798347 0.0175438597798347 0.0175438597798347;0.0233918130397797 0.0233918130397797 0.0233918130397797;0.0292397662997246 0.0292397662997246 0.0292397662997246;0.0350877195596695 0.0350877195596695 0.0350877195596695;0.0409356728196144 0.0409356728196144 0.0409356728196144;0.0467836260795593 0.0467836260795593 0.0467836260795593;0.0526315793395042 0.0526315793395042 0.0526315793395042;0.0584795325994492 0.0584795325994492 0.0584795325994492;0.0643274858593941 0.0643274858593941 0.0643274858593941;0.070175439119339 0.070175439119339 0.070175439119339;0.0760233923792839 0.0760233923792839 0.0760233923792839;0.0818713456392288 0.0818713456392288 0.0818713456392288;0.0877192988991737 0.0877192988991737 0.0877192988991737;0.0935672521591187 0.0935672521591187 0.0935672521591187;0.0994152054190636 0.0994152054190636 0.0994152054190636;0.105263158679008 0.105263158679008 0.105263158679008;0.111111111938953 0.111111111938953 0.111111111938953;0.116959065198898 0.116959065198898 0.116959065198898;0.122807018458843 0.122807018458843 0.122807018458843;0.128654971718788 0.128654971718788 0.128654971718788;0.134502917528152 0.134502917528152 0.134502917528152;0.140350878238678 0.140350878238678 0.140350878238678;0.146198838949203 0.146198838949203 0.146198838949203;0.152046784758568 0.152046784758568 0.152046784758568;0.157894730567932 0.157894730567932 0.157894730567932;0.163742691278458 0.163742691278458 0.163742691278458;0.169590651988983 0.169590651988983 0.169590651988983;0.175438597798347 0.175438597798347 0.175438597798347;0.181286543607712 0.181286543607712 0.181286543607712;0.187134504318237 0.187134504318237 0.187134504318237;0.192982465028763 0.192982465028763 0.192982465028763;0.198830410838127 0.198830410838127 0.198830410838127;0.204678356647491 0.204678356647491 0.204678356647491;0.210526317358017 0.210526317358017 0.210526317358017;0.216374278068542 0.216374278068542 0.216374278068542;0.222222223877907 0.222222223877907 0.222222223877907;0.253333330154419 0.253333330154419 0.253333330154419;0.284444451332092 0.284444451332092 0.284444451332092;0.315555542707443 0.315555542707443 0.315555542707443;0.346666663885117 0.346666663885117 0.346666663885117;0.37777778506279 0.37777778506279 0.37777778506279;0.408888876438141 0.408888876438141 0.408888876438141;0.439999997615814 0.439999997615814 0.439999997615814;0.471111118793488 0.471111118793488 0.471111118793488;0.502222239971161 0.502222239971161 0.502222239971161;0.533333361148834 0.533333361148834 0.533333361148834;0.564444422721863 0.564444422721863 0.564444422721863;0.595555543899536 0.595555543899536 0.595555543899536;0.626666665077209 0.626666665077209 0.626666665077209;0.657777786254883 0.657777786254883 0.657777786254883;0.688888907432556 0.688888907432556 0.688888907432556;0.720000028610229 0.720000028610229 0.720000028610229;0.751111090183258 0.751111090183258 0.751111090183258;0.782222211360931 0.782222211360931 0.782222211360931;0.813333332538605 0.813333332538605 0.813333332538605;0.844444453716278 0.844444453716278 0.844444453716278;0.875555574893951 0.875555574893951 0.875555574893951;0.906666696071625 0.906666696071625 0.906666696071625;0.937777757644653 0.937777757644653 0.937777757644653;0.968888878822327 0.968888878822327 0.968888878822327;1 1 1]);

% Create axes
axes1 = axes('Parent',figure1,'Position',[0.13 0.11 0.4090625 0.815],...
    'FontSize',14);
% Uncomment the following line to preserve the X-limits of the axes
 xlim(axes1,[0.7 1.5]);
% Uncomment the following line to preserve the Y-limits of the axes
 ylim(axes1,[0 0.9375]);
hold(axes1,'all');

% Create contour
contour(xdata1,ydata1,zdata1,'LineColor',[0 0 0],'Fill','on','Parent',axes1);

% Create quiver
quiver(xdata1,ydata1,U1,V1,'LineWidth',2,...
    'AutoScaleFactor',0.189999997615814,...
    'Parent',axes1);

% Create xlabel
xlabel('R (m)','FontSize',20);

% Create ylabel
ylabel('Z (m)','FontSize',20);

% Create zlabel
zlabel({'',''},'Visible','off');

% Create colorbar
colorbar('peer',axes1,'FontSize',14);
% Resize the axes in order to prevent it from shrinking.
set(axes1,'Position',[0.13 0.11 0.4090625 0.815]);
