%clear;
%clc;

figure
subplot(2,1,1);
load('Fig8_Case1.mat');
POT=[POT1;POT2;POT3;POT4];
IndexStable = find(POT(:,4)>0);
IndexUnstable = find(POT(:,4)<0);
scatter(POT(IndexStable,1),POT(IndexStable,2), 'o','blue','filled','SizeData',25);
hold on
scatter(POT(IndexUnstable,1),POT(IndexUnstable,2), 'o','magenta','filled','SizeData',25);
%xlabel('\omega_B (admin)');
%xlim([1.25 1.8])
ylabel('H (admin)');
title('Case #1')
box('on');
%grid on
%set(gca,'FontSize',20,'FontWeight','bold','LineWidth',2)

subplot(2,1,2); 
load('Fig8_Case2.mat');
POT = [POT1;POT2;POT3;POT4;POT5];
IndexStable = find(POT(:,4)>0);
IndexUnstable = find(POT(:,4)<0);
scatter(POT(IndexStable,1),POT(IndexStable,2), 'o','blue','filled','SizeData',25);
hold on
scatter(POT(IndexUnstable,1),POT(IndexUnstable,2), 'o','magenta','filled','SizeData',25);
xlabel('\omega_B (admin)');
%xlim([1.25 1.8])
ylabel('H (admin)');
title('Case #2')
box('on');
%grid on
%set(gca,'FontSize',20,'FontWeight','bold','LineWidth',2)