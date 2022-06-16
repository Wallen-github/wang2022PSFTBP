%clear;
%clc;
YHC4 = load('YHC4.DAT');

figure
scatter(YHC4(:,9),YHC4(:,11), '.','black');
xlabel('t (adim)');
ylabel('Momentum (adim)');
box('on');
grid on
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)

figure
scatter(YHC4(:,9),YHC4(:,10), '.','black');
xlabel('t (adim)');
ylabel('Energy (adim)');
box('on');
grid on
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)


figure;
subplot(2,2,1);
plot(YHC4(:,1).*cos(YHC4(:,2)),YHC4(:,1).*sin(YHC4(:,2)),'color','black','LineWidth',1)
%title('\omega_B=1.4956');
grid on
xlabel('x (adim)');
ylabel('y (adim)');

subplot(2,2,2);
plot(YHC4(:,7)+YHC4(:,8),YHC4(:,1),'color','black','LineWidth',1);
grid on
xlabel('\theta_B1 (adim)');
ylabel('S (adim)');

subplot(2,2,3);
plot(YHC4(:,9),mod(YHC4(:,2),2*pi),'color','black','LineWidth',1);
grid on
ylabel('\theta (deg)');
xlabel('t (adim)');

subplot(2,2,4);
plot(YHC4(:,9),mod(YHC4(:,4),pi*2),'color','black','LineWidth',1);
grid on
xlabel('t (adim)');
ylabel('\theta_A (\times 2\pi)');

figure
scatter(YHC4(:,7)+YHC4(:,8),YHC4(:,2), '.','black');
xlabel('\omega_B');
ylabel('\theta');
box('on');
grid on
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)

figure
plot(YHC4(:,7)+YHC4(:,8),YHC4(:,1),'color','black','LineWidth',1);
xlabel('\omega_B (adim)');
ylabel('S (adim)');
box('on');
grid on
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',2)

