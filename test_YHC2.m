close all;
%clear all;
YHC2 = load('YHC2.DAT');

figure
plot(YHC2(:,9),YHC2(:,10), 'LineWidth', 2, 'Color', 'black','DisplayName','3rd-body');
xlabel('t');
ylabel('Energy');
box('on');
grid on
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)


figure
scatter(YHC2(:,2),YHC2(:,3), '.','black');
ylabel('delta1');
xlabel('phi');
box('on');
grid on
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)

figure
scatter(YHC2(:,5).*cos(YHC2(:,6)),YHC2(:,5).*sin(YHC2(:,6)), '.','black');
ylabel('x');
xlabel('y');
box('on');
grid on
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)