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

% Curve fit all data:
x = [x1,x2,x3,x4,x5,x6,x7,x8,x9];
y = [y1,y2,y3,y4,y5,y6,y7,y8,y9];
IC = [70,2];
plotData(x,y,IC,'All talkers, all rooms');

% Curve fit anechoic only:
x = [x1,x2,x3];
y = [y1,y2,y3];
IC = [70,2];
plotData(x,y,IC,'All talkers, no RV');

% Curve fit room 10:
x = [x4,x5,x6];
y = [y4,y5,y6];
IC = [70,2];
plotData(x,y,IC,'All talkers, low RV');

% Curve fit reverb room:
x = [x7,x8,x9];
y = [y7,y8,y9];
IC = [70,2];
plotData(x,y,IC,'All talkers, high RV');

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

% Curve fit quiet, low BGN & low RV, high BGN & high RV:
x = [x1,x5,x9];
y = [y1,y5,y9];
IC = [70,2];
plotData(x,y,IC,'All talkers, NoRV&NoBGN, LowRV&LowBGN, HighRV&HighBGN');



return;

function plotData(x,y,IC,titleText)

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
set([gca,get(gca,'xlabel'),get(gca,'ylabel'),get(gca,'title')],'fontname','arial','fontsize',12);

return;

% Bye!