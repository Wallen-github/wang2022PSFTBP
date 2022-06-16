%clear;
%clc;

%POT = load('PO.DAT');
%POT = POT_Case1_left2;

IndexStable = find(POT(:,4)>0);
IndexUnstable = find(POT(:,4)<0);

figure
scatter(POT(IndexStable,1),POT(IndexStable,2), '.','blue');
hold on
scatter(POT(IndexUnstable,1),POT(IndexUnstable,2), '.','magenta');
xlabel('\omega_B (admin)');
ylabel('H (admin)');
box('on');
grid on
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)

figure
scatter(POT(IndexStable,1),POT(IndexStable,5), '.','blue');
hold on
scatter(POT(IndexUnstable,1),POT(IndexUnstable,5), '.','magenta');
xlabel('\omega_B (admin)');
ylabel('S (admin)');
box('on');
grid on
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)

figure
scatter(POT(IndexStable,1),POT(IndexStable,3), '.','blue');
hold on
scatter(POT(IndexUnstable,1),POT(IndexUnstable,3), '.','magenta');
xlabel('\omega_B (admin)');
ylabel('Energy');
box('on');
grid on
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)


