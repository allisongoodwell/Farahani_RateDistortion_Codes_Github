%figure 3 : 5/5/21
% Illustrate original Ta forcing data and simplified Ta 
% to N=2, 3, 4, and 5 levels  of quantization for growing season.

clear 
close all

load T
load Ta_crop
load quantdata

steps_per_day = 4*24; %15 minute data, each day
step = 1/steps_per_day; %15 minute step
X=0.01041:step:365;
figure(1)
subplot(4,2,[1,3,5,7]);
plot(X,Ta_crop,'b','LineWidth',1)
hold on
xline(X(1,14305),'r','LineWidth',1);
xline(X(1,24000),'r','LineWidth',1);
annotation('rectangle', [0.266, 0.62, 0.0935, 0.22],'Color','r')
xlabel('Day of Year', 'FontSize',12,'FontWeight','bold') 
ylabel({'a)';'Temperature (\circ C) '}, 'FontSize',12,'FontWeight','bold')
xlim([1 365]);
hold off
 
subplot(4,2,2);
plot(X(1,14305:24000),Ta_Q(:,1),'b','LineWidth',1)
ax = gca;
ax.YColor = 'r';
ax.YAxis.LineWidth = 1.5;
ylabel({'b)';'N=2'}, 'FontSize',12,'FontWeight','bold','Color','k')
ylim([15 30])
xlim([150 250]);

subplot(4,2,4);
plot(X(1,14305:24000),Ta_Q(:,2),'b','LineWidth',1)
ax = gca;
ax.YColor = 'r';
ax.YAxis.LineWidth = 1.5;
ylabel({'c)';'N=3'}, 'FontSize',12,'FontWeight','bold','Color','k')
ylim([15 30])
xlim([150 250]);

subplot(4,2,6);
plot(X(1,14305:24000),Ta_Q(:,3),'b','LineWidth',1)
ax = gca;
ax.YColor = 'r';
ax.YAxis.LineWidth = 1.5;
ylabel({'d)';'N=4'}, 'FontSize',12,'FontWeight','bold','Color','k')
ylim([15 30])
xlim([150 250]);

subplot(4,2,8);
plot(X(1,14305:24000),Ta_Q(:,4),'b','LineWidth',1)
xlabel('Day of Year', 'FontSize',12,'FontWeight','bold')
ax = gca;
ax.YColor = 'r';
ax.YAxis.LineWidth = 1.5;
ylabel({'e)';'N=5'}, 'FontSize',12,'FontWeight','bold','Color','k')
ylim([15 30])
xlim([150 250]);