%% Figure setting
format compact;
set(0, "defaultFigureWindowStyle", "normal");

%% font
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');


% GUI font
set(0, 'defaultUicontrolFontName', 'Meiryo UI');
% Axis font
set(0, 'defaultAxesFontName', 'Times');
% Title, legend font
set(0, 'defaultTextFontName', 'Times');

% % GUI fontsize
% set(0, 'defaultUicontrolFontSize', 9);
% 
% Axis fontsize
set(0, 'defaultAxesFontSize', 15);
set(0, 'defaultLegendFontSize', 15);
set(0, 'defaultColorbarFontSize', 15);
set(0, 'defaultTextarrowshapeFontSize', 8);
set(0, 'defaultTextboxshapeFontSize', 8);


% Title, legend fontsize
set(0, 'defaultTextFontSize', 10);

%% other
% axis line thickness
set(0, 'DefaultAxesLineWidth', 2);
set(0, 'DefaultLegendAutoUpdate', 'off');

% legend object line thickness
set(0, 'DefaultLineLineWidth', .01);
set(0, 'DefaultLineMarkerSize', 10);
set(0, 'DefaultStairLineWidth', 2);
set(0, 'DefaultStairMarkerSize', 10);

% default color map
cmap = spring(128);
set(0, 'defaultFigureColormap', cmap);

% plot color default
corder = [0,0,1;1,0,0;0.8,0.8,0;0.66015625,0.66015625,0.66015625;0,0,0;1,0.64453125,0;1,0,1;0,0.5,0.5;0,0,0.54296875;0,0.390625,0;0,1,1;0.59765625,0.1953125,0.796875];
set(0, 'defaultAxesColorOrder', corder);

% figure size
set(0, 'defaultFigureUnits','pixels')
set(0, 'defaultFigurePosition',[100 100 650 400])

clearvars cmap corder;