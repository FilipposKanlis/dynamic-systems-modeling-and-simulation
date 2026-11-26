% Filippos Kanlis 9829
% Ergasia 2
% thema 3

clear;
clc;
close all;

a11 = -0.25;
a12 = -3;
a21 = 5;
a22 = 0;
b1 = 0.5;
b2 = 1.5;
theta_m = 0.01;

u = @(t) 3.5*sin(7.2*t) + 2*sin(11.7*t);
t= 0:0.01:20;

% Estimate g1 and g2
select_method=1;
[g1_best,g2_best] = estimate_g1_g2(select_method,a11,a12,a21,a22,b1,b2,theta_m,u,t);

% Eksiswseis katastashs, epilush me g1_best, g2_best
[t,x]=solveDiffEq(a11,a12,a21,a22,b1,b2,theta_m,u,t,g1_best,g2_best);


% Απεικόνιση x1,x2 και εκτιμήσεων x1,x2.
fig=figure();
subplot(2,1,1)
hold on;
plot(t,x(:,1));
plot(t,x(:,9));
title('[Mixed structure] $x_1$ and $\hat{x_1}$','interpreter','latex','FontSize',16);
xlabel('Time [s]','FontSize',15);
legend('$x_1$','$\hat{x_1}$','interpreter','latex','FontSize',10);

subplot(2,1,2)
hold on;
plot(t,x(:,2));
plot(t,x(:,10));
title('[Mixed structure] $x_2$ and $\hat{x_2}$','interpreter','latex','FontSize',16);
xlabel('Time [s]','FontSize',15);
legend('$x_2$','$\hat{x_2}$','interpreter','latex','FontSize',10);
figname="x_1&x_2.jpg";
saveas(fig,figname)

% Απεικόνιση σφάλματων e1 και e2, όπου ei = xiReal - xiModel με i=1,2.
fig=figure();
subplot(2,1,1)
plot(t,x(:,1)- x(:,9));
title('[Mixed structure] $e_1$ = $x_1$ - $\hat{x_1}$','interpreter','latex','FontSize',16);
xlabel('Time [s]','FontSize',15);

subplot(2,1,2)
plot(t,x(:,2)- x(:,10));
title('[Mixed structure] $e_2$ = $x_2$ - $\hat{x_2}$','interpreter','latex','FontSize',16);
xlabel('Time [s]','FontSize',15);
figname="e_1&e_2.jpg";
saveas(fig,figname)

% Απεικόνιση παραμέτρων a11,a12,a21,a22 και των εκτιμήσεων τους.
fig=figure();
subplot(2,2,1)
hold on;
plot(t,x(:,3));
yline(a11,'-k');
hold off
title('$a_{11}$ and $\hat{a_{11}}$','interpreter','latex','FontSize',16);
legend('$\hat{a_{11}}$','real $a_{11}$','interpreter','latex','FontSize',10);
xlabel('Time [s]','FontSize',15);

subplot(2,2,2)
hold on;
plot(t,x(:,4));
yline(a12,'-k');
hold off
title('$a_{12}$ and $\hat{a_{12}}$','interpreter','latex','FontSize',16);
legend('$\hat{a_{12}}$','real $a_{12}$','interpreter','latex','FontSize',10);
xlabel('Time [s]','FontSize',15);

subplot(2,2,3)
hold on;
plot(t,x(:,5));
yline(a21,'-k');
hold off
title('$a_{21}$ and $\hat{a_{21}}$','interpreter','latex','FontSize',16);
legend('$\hat{a_{21}}$','real $a_{21}$','interpreter','latex','FontSize',10);
xlabel('Time [s]','FontSize',15);

subplot(2,2,4)
hold on;
plot(t,x(:,6));
yline(a22,'-k');
hold off
title('$a_{22}$ and $\hat{a_{22}}$','interpreter','latex','FontSize',16);
legend('$\hat{a_{22}}$','real $a_{22}$','interpreter','latex','FontSize',10);
xlabel('Time [s]','FontSize',15);
figname="a_xx.jpg";
saveas(fig,figname)

% Απεικόνιση παραμέτρων b1,b2 και των εκτιμήσεων τους.
fig=figure();
subplot(2,1,1)
hold on;
plot(t,x(:,7));
yline(b1,'-k');
hold off
title('$b_{1}$ and $\hat{b_{1}}$','interpreter','latex','FontSize',16);
legend('$\hat{b_{1}}$','real $b_{1}$','interpreter','latex','FontSize',10);
xlabel('Time [s]','FontSize',15);

subplot(2,1,2)
hold on;
plot(t,x(:,8));
yline(b2,'-k');
hold off
title('$b_{2}$ and $\hat{b_{2}}$','interpreter','latex','FontSize',16);
legend('$\hat{b_{2}}$','real $b_{2}$','interpreter','latex','FontSize',10);
xlabel('Time [s]','FontSize',15);
figname="b1&b2.jpg";
saveas(fig,figname)
 

function [t,x]=solveDiffEq(a11,a12,a21,a22,b1,b2,theta_m,u,t,gi,gj)
                odefun = @(t,x) [a11*x(1) + a12*x(2) + b1*u(t);
                                a21*x(1) + a22*x(2) + b2*u(t);
                                gi*x(9)*(x(1)-x(9));
                                gi*x(10)*(x(1)-x(9));
                                gi*x(9)*(x(2)-x(10));
                                gi*x(10)*(x(2)-x(10));
                                gj*u(t)*(x(1)-x(9));
                                gj*u(t)*(x(2)-x(10));
                                x(3)*x(9) + x(4)*x(10) + x(7)*u(t) + theta_m*(x(1)-x(9));
                                x(5)*x(9) + x(6)*x(10) + x(8)*u(t) + theta_m*(x(2)-x(10))];
                [t,x] = ode45(odefun,t,[0,0,0,0,0,0,0,0,0,0]);
end

function [g1_best,g2_best] = estimate_g1_g2(select_method,a11,a12,a21,a22,b1,b2,theta_m,u,t)

    % select_method = 1 : Ektimhsh twn g1 kai g2 gia parallhlh domh
    % select_method = 2: Ektimhsh twn g1 kai g2 gia mikth domh
    
    if select_method == 1
        fprintf("[Mixed Structure] Find optimal g1 and g2\n");
        g1_range=20:5:50;
        g2_range=20:5:50; %autes oi times exoun stathera xaraxthkristika: grhgorh statheropoihsh, grhgorh aposvesh talantwshs 
        min = 1;
        for gj = g2_range
            for gi = g1_range
                [t,x]=solveDiffEq(a11,a12,a21,a22,b1,b2,theta_m,u,t,gi,gj);
                mas = (sum(abs(x(:,1)- x(:,9))))/length(x); % mas = meso apoluto sfalma 
                if mas < min
                    min = mas;
                    g1_best = gi;
                    g2_best = gj;
                end
                fprintf("for g1 = %d and g2 = %d Meso Apoluto sfalma = %f\n",gi,gj,mas)
            end
        end
        fprintf("Optimal pair: [g1,g2] = [%d,%d]\n",g1_best,g2_best);
    end
end