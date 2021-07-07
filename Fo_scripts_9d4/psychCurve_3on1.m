function psychCurve()
% This function is amazing.
% psychCurve
% Lisa Kopf, 2/11/2014
% Reference: http://courses.washington.edu/matlab1/Lesson_5.html

close all;

% Define x-values & y-values (BGN intensity (x), intelligibility (y), anechoic chamber)
%No BGN
x1 = [0 0 0 0 0 0];
y1 = [0	0.158730159	0.584795322	0	5.333333333	2.189542484];
%Low BGN
x2 = [40    40  40  40  40  40];
y2 = [0	0	0.303030303	2.083333333	4.460784314	3.518518519];
% High BGN
x3 = [50    60  60  60  60  60];
y3 = [1.388888889	14.95726496	39.64506173	1.666666667	79.93464052	65];

% Define x-values & y-values (BGN intensity (x), intelligibility (y), room 10)
%No BGN
x4 = [0 0 0 0 0 0];
y4 = [0	1.944444444	0.245098039	2.777777778	6.888888889	2.407407407];
%Low BGN
x5 = [40    40  40  40  40  40];
y5 = [0.681818182	3.260869565	3.58974359	4.393939394	8.148148148	25.80246914];
% High BGN
x6 = [50    60  60  60  60  60];
y6 = [1.481481481	48.47953216	20.67460317	19.29292929	79.6969697	69.73958333];

% Define x-values & y-values (BGN intensity (x), intelligibility (y), reverb room)
%No BGN
x7 = [0 0 0 0 0 0];
y7 = [2.777777778	4.074074074	5.802469136	12.28395062	19.20289855	22.77777778];
%Low BGN
x8 = [40    40  40  40  40  40];
y8 = [6.068376068	28.125	5.333333333	5.882352941	30.36111111	25.55555556];
% High BGN
x9 = [50    50  60  60  60  60];
y9 = [4.965277778	16.94444444	79.86111111	47.32638889	95	80.41666667];


% Define x-values & y-values (RV time (x), intelligibility (y), no BGN)
%No RV
x11 = [0 0 0 0 0 0];
y11 = [0	0.158730159	0.584795322	0	5.333333333	2.189542484];
%Low RV
x12 = [0.8  0.8 0.8 0.8 0.8 0.8];
y12 = [0	1.944444444	0.245098039	2.777777778	6.888888889	2.407407407];
% High RV
x13 = [2.5  2.5 2.5 2.5 2.5 2.5];
y13 = [2.777777778	4.074074074	5.802469136	12.28395062	19.20289855	22.77777778];

% Define x-values & y-values (RV time (x), intelligibility (y), low BGN)
%No RV
x14 = [0 0 0 0 0 0];
y14 = [0.00	0.00	0.30	2.08	4.46	3.52];
%Low RV
x15 = [0.8  0.8 0.8 0.8 0.8 0.8];
y15 = [0.68	3.26	3.59	4.39	8.15	25.80];
% High RV
x16 = [2.5  2.5 2.5 2.5 2.5 2.5];
y16 = [6.07	28.13	5.33	5.88	30.36	25.56];

% Define x-values & y-values (RV time (x), intelligibility (y), high BGN)
%No RV
x17 = [0 0 0 0 0 0];
y17 = [1.39	14.96	39.65	1.67	79.93	65.00];
%Low RV
x18 = [0.8  0.8 0.8 0.8 0.8 0.8];
y18 = [1.48	48.48	20.67	19.29	79.70	69.74];
% High RV
x19 = [2.5  2.5 2.5 2.5 2.5 2.5];
y19 = [4.97	16.94	79.86	47.33	95.00	80.42];

% % Curve fit all data:
% x = [x1,x2,x3,x4,x5,x6,x7,x8,x9];
% y = [y1,y2,y3,y4,y5,y6,y7,y8,y9];
% IC = [70,2];
% plotData(x,y,IC,'All talkers, all rooms');

% Curve fit anechoic only:
x = [x1,x2,x3];
y = [y1,y2,y3];
IC = [70,2];
outputAll = logisticCurveFit(x,y,IC);
figure;
hold on;
plot(x,y,'bo','MarkerFaceColor','b');
xModel = [20:100];
yModel = feval(outputAll.fitResults,xModel);
plot(xModel,yModel,'b-','LineWidth',2);
xlabel('Noise level, dB');
ylabel('Percent error');
ylim([0 100]);
grid on;
set([gca,get(gca,'xlabel'),get(gca,'ylabel')],'fontname','arial','fontsize',12);

hold on;

% Curve fit room 10:
x = [x4,x5,x6];
y = [y4,y5,y6];
IC = [70,2];
outputSome = logisticCurveFit(x,y,IC);
plot(x,y,'go','MarkerFaceColor','g');
xModel = [20:100];
yModel = feval(outputSome.fitResults,xModel);
plot(xModel,yModel,'g-','LineWidth',2);
xlabel('Noise level, dB');
ylabel('Percent error');
ylim([0 100]);
grid on;
set([gca,get(gca,'xlabel'),get(gca,'ylabel')],'fontname','arial','fontsize',12);

% Curve fit reverb room:
x = [x7,x8,x9];
y = [y7,y8,y9];
IC = [70,2];
outputNone= logisticCurveFit(x,y,IC);
plot(x,y,'ro','MarkerFaceColor','r');
xModel = [20:100];
yModel = feval(outputNone.fitResults,xModel);
plot(xModel,yModel,'r-','LineWidth',2);
xlabel('Noise level, dB');
ylabel('Percent error');
ylim([0 100]);
grid on;
set([gca,get(gca,'xlabel'),get(gca,'ylabel')],'fontname','arial','fontsize',12);


% Curve fit no BGN:
x = [x11,x12,x13];
y = [y11,y12,y13];
IC = [70,2];
outputAll = logisticCurveFit(x,y,IC);
figure;
hold on;
plot(x,y,'bo','MarkerFaceColor','b');
xModel = [-10:10];
yModel = feval(outputAll.fitResults,xModel);
plot(xModel,yModel,'b-','LineWidth',2);
xlabel('Reverberation Time, sec');
ylabel('Percent error');
ylim([0 100]);
grid on;
set([gca,get(gca,'xlabel'),get(gca,'ylabel')],'fontname','arial','fontsize',12);

hold on;

% Curve fit low BGN:
x = [x14,x15,x16];
y = [y14,y15,y16];
IC = [70,2];
outputSome = logisticCurveFit(x,y,IC);
plot(x,y,'go','MarkerFaceColor','g');
xModel = [-10:10];
yModel = feval(outputSome.fitResults,xModel);
plot(xModel,yModel,'g-','LineWidth',2);
ylim([0 100]);
grid on;
set([gca,get(gca,'xlabel'),get(gca,'ylabel')],'fontname','arial','fontsize',12);

% Curve fit high BGN:
x = [x17,x18,x19];
y = [y17,y18,y19];
IC = [70,2];
outputNone= logisticCurveFit(x,y,IC);
plot(x,y,'ro','MarkerFaceColor','r');
xModel = [-10:10];
yModel = feval(outputNone.fitResults,xModel);
plot(xModel,yModel,'r-','LineWidth',2);
ylim([0 100]);
grid on;
set([gca,get(gca,'xlabel'),get(gca,'ylabel')],'fontname','arial','fontsize',12);

%  % Curve fit quiet:
% x = [x1,x4,x7];
% y = [y1,y4,y7];
% IC = [70,2];
% plotData(x,y,IC,'All talkers, no BGN');
% 
% % Curve fit low BGN:
% x = [x2,x5,x8];
% y = [y2,y5,y8];
% IC = [70,2];
% plotData(x,y,IC,'All talkers, low BGN');
% 
% % Curve fit high BGN:
% x = [x3,x6,x9];
% y = [y3,y6,y9];
% IC = [70,2];
% plotData(x,y,IC,'All talkers, high BGN');

% % Curve fit quiet, low BGN & low RV, high BGN & high RV:
x = [x1,x5,x9];
y = [y1,y5,y9];
IC = [70,2];
outputAll = logisticCurveFit(x,y,IC);
figure;
hold on;
plot(x,y,'ko');
xModel = [20:100];
yModel = feval(outputAll.fitResults,xModel);
plot(xModel,yModel,'r-','LineWidth',2);
xlabel('Noise level, dB');
ylabel('Percent error');
ylim([0 100]);
grid on;
title(titleText);
set([gca,get(gca,'xlabel'),get(gca,'ylabel')],'fontname','arial','fontsize',12);



return;

% function plotData(x,y,IC,titleText)
% 
% outputAll = logisticCurveFit(x,y,IC);
% figure;
% hold on;
% plot(x,y,'ko');
% xModel = [20:100];
% yModel = feval(outputAll.fitResults,xModel);
% plot(xModel,yModel,'r-','LineWidth',2);
% xlabel('Noise level, dB');
% ylabel('Percent error');
% ylim([0 100]);
% grid on;
% title(titleText);
% set([gca,get(gca,'xlabel'),get(gca,'ylabel'),get(gca,'title')],'fontname','arial','fontsize',12);
% 
% return;

% Bye!