% Filippos Kanlis 9829
% Ergasia 2 
% thema 2

clear;
clc;
close all;

a = 3;
b = 0.5;
h0 = 0.5;
f = 40;
theta_m = 0.1;

u = @(t) 10*sin(3*t);
h = @(t,f,h0)h0*sin(2*pi*f*t);
t = 0:0.01:41.89;

% select_method = 1 : Ektimhsh twn g1 kai g2 gia parallhlh domh
% select_method = 2: Ektimhsh twn g1 kai g2 gia mikth domh

% Estimate g1 and g2
select_method=1;
[g1_best,g2_best] = estimate_g1_g2(select_method,a,b,theta_m,u,t);

%Parameter estimtion with f = 40 Hz and h0 = 0.5
%Ektimisi parametrwn me xrhsh parallhlhs domhs
fi=[10,20,30,40,50];
h0i=[0.5,0.75,1,1.5];

for f=fi
    h0=0.5;
    [x_pred,a_pred,b_pred] = parallel_structure(a,b,g1_best,g2_best,u,h,f,h0,t);
end

for h0=h0i
    f=40;
    [x_pred,a_pred,b_pred] = parallel_structure(a,b,g1_best,g2_best,u,h,f,h0,t);
end

% Estimate g1 and g2
select_method=2;
[g1_best,g2_best] = estimate_g1_g2(select_method,a,b,theta_m,u,t);

%Parameter estimtion with f = 40 Hz and h0 = 0.5
%Ektimisi parametrwn me xrhsh mikths domhs
fi=[10,20,30,40,50];
h0i=[0.5,0.75,1,1.05];

for f=fi
    h0=0.5;
    [x_pred,a_pred,b_pred] = mixed_structure(a,b,theta_m,g1_best,g2_best,u,h,f,h0,t);
end

for h0=h0i
    f=40;
    [x_pred,a_pred,b_pred] = mixed_structure(a,b,theta_m,g1_best,g2_best,u,h,f,h0,t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% functions   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [g1_best,g2_best] = estimate_g1_g2(select_method,a,b,theta_m,u,t)
    
    % select_method = 1 : Ektimhsh twn g1 kai g2 gia parallhlh domh
    % select_method = 2: Ektimhsh twn g1 kai g2 gia mikth domh
    
    if select_method == 1
        fprintf("[Parallel Structure] Find optimal g1 and g2\n");
        g1_range=5:10:30;
        g2_range=[5,20,50]; %autes oi times exoun stathera xaraxthkristika: grhgorh statheropoihsh, grhgorh aposvesh talantwshs 
        min = 1;
        for gj = g2_range
            for gi = g1_range
                odefun = @(t,x) [-a*x(1)+ b*u(t);
                -gi*(x(1)-x(4))*x(4);
                 gj*(x(1)-x(4))*u(t);
                -x(2)*x(4)+x(3)*u(t)];
                [t,x] = ode45(odefun,t,[0,0,0,0]);
                mas = (sum(abs(x(:,1)- x(:,4))))/length(x); % mas = meso apoluto sfalma 
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
    
    if select_method == 2
        fprintf("[Mixed Structure] Find optimal g1 and g2\n");
        g1_range=5:5:15;
        g2_range=10:5:30;
        min = 1;
        for gj = g2_range
            for gi = g1_range
                 odefun = @(t,x) [-a*x(1)+ b*u(t);
                -gi*(x(1)-x(4))*x(4);
                gj*(x(1)-x(4))*u(t);
                -x(2)*x(4)+x(3)*u(t)+theta_m*(x(1)-x(4))];
                [t,x] = ode45(odefun,t,[0,0,0,0]);
                mas = (sum(abs(x(:,1)- x(:,4))))/length(x); % mas = meso apoluto sfalma 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% parallel_structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_pred,a_pred,b_pred] = parallel_structure(a,b,g1_best,g2_best,u,h,f,h0,t)    
    %xwris thorivo
    odefun = @(t,x) [-a*x(1)+ b*u(t);
                     -g1_best*(x(1)-x(4))*x(4);
                      g2_best*(x(1)-x(4))*u(t);
                      -x(2)*x(4)+x(3)*u(t)];
    [t,x] = ode45(odefun,t,[0,0,0,0]);
    %me thorivo
    odefun_noise = @(t,x) [-a*x(1)+ b*u(t);
                     -g1_best*(x(1)+h(t,f,h0)-x(4))*x(4);
                      g2_best*(x(1)+h(t,f,h0)-x(4))*u(t);
                      -x(2)*x(4)+x(3)*u(t)];
    [t,xnoise] = ode45(odefun_noise,t,[0,0,0,0]);

    x_pred = xnoise(:,4);
    a_pred = xnoise(:,2);
    b_pred = xnoise(:,3);
    mas = (sum(abs(x(:,1)- x(:,4))))/length(x); 
    mas_noise = (sum(abs(xnoise(:,1)- xnoise(:,4))))/length(x);

    %x and x_pred with and without noise
    fig=figure();
    subplot(2,1,1)
    hold on;
    plot(t,x(:,1));
    plot(t,x(:,4));
    hold off
    grid on;
    xticks(0:5:60)
    title('[Parallel structure] $x$ and $\hat{x}$ without noise','interpreter','latex','FontSize',16);
    xlabel('Time [s]','FontSize',15);
    legend('$x$','$\hat{x}$','interpreter','latex','FontSize',16);
    
    subplot(2,1,2)
    hold on;
    plot(t,xnoise(:,1));
    plot(t,xnoise(:,4));
    hold off
    grid on;
    xticks(0:5:60)
    title({'[Parallel structure] $x$ and $\hat{x}$ with noise';['f =',num2str(f),', h0 =',num2str(h0)]},'interpreter','latex','FontSize',16);
    xlabel('Time [s]','FontSize',15);
    legend('$x$','$\hat{x}$','interpreter','latex','FontSize',16);
    figname="ps f="+f+ "_h0="+ h0+"_x&xhat.jpg";
    saveas(fig,figname)
    
                    
    % e = x - x_pred with and without noise
    fig=figure();
    subplot(2,1,1)
    plot(t,x(:,1)- x(:,4));
    title('[Parallel structure] e = $x$ - $\hat{x}$ without noise','interpreter','latex','FontSize',16);
    xlabel('Time [s]','FontSize',15);
    grid on;
    xticks(0:5:60)
    dim = [.55 .6 .35 .07];
    str = strcat('Mean Absolute Error =',num2str(mas));
    annotation('textbox',dim,'String',str,'FontSize',10)
    
    subplot(2,1,2)
    plot(t,xnoise(:,1)- xnoise(:,4));
    title({'[Parallel structure] e = $x$ - $\hat{x}$ with noise';['f =',num2str(f),', h0 =',num2str(h0)]},'interpreter','latex','FontSize',16);
    xlabel('Time [s]','FontSize',15);
    grid on;
    xticks(0:5:60)
    dim = [.55 .125 .35 .07];
    str = strcat('Mean Absolute Error =',num2str(mas_noise));
    annotation('textbox',dim,'String',str,'FontSize',10)
    figname="ps f="+f+ "_h0="+ h0+"_x-xhat.jpg";
    saveas(fig,figname)
   
    
    
    % a and b with and without noise
    fig=figure();
    subplot(2,1,1)
    hold on;
    plot(t,xnoise(:,2));
    plot(t,x(:,2));
    yline(1.5,'-k');
    hold off
    title({'[Parallel structure] a and $\hat{a}$ with and without noise';['f =',num2str(f),', h0 =',num2str(h0)]},'interpreter','latex','FontSize',16);
    legend('$\hat{a}$ with noise','$\hat{a}$ without noise','$a_{real}$','interpreter','latex');
    xlabel('Time [s]','FontSize',15);
    grid on;
    xticks(0:5:60)
    
    subplot(2,1,2)
    hold on;
    plot(t,xnoise(:,3));
    plot(t,x(:,3));
    yline(2,'-k');
    hold off
    title({'[Parallel structure] b and $\hat{b}$ with and without noise';['f =',num2str(f),', h0 =',num2str(h0)]},'interpreter','latex','FontSize',16);
    legend('$\hat{b}$ with noise','$\hat{b}$ without noise','$b_{real}$','interpreter','latex');
    xlabel('Time [s]','FontSize',15);
    grid on;
    xticks(0:5:60)
    figname="ps f="+f+ "_h0="+ h0+"_a&b.jpg";
    saveas(fig,figname)
   
    fprintf("[Parallel Structure] f = %d and h = %d: Mean Absolute Error = %f\n",f,h0,mas_noise);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% mixed_structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_pred,a_pred,b_pred] = mixed_structure(a,b,theta_m,g1_best,g2_best,u,h,f,h0,t)
    odefun = @(t,x) [-a*x(1) + b*u(t);
                         -g1_best*(x(1)-x(4))*x(1);
                          g2_best*(x(1)-x(4))*u(t);
                          -x(2)*x(4)+x(3)*u(t)+theta_m*(x(1)-x(4))];
    
    % ektimhsh twn katastasewn tou susthmatos otan exoume thorivo
    odefun_noise = @(t,x) [-a*x(1) + b*u(t);
                         -g1_best*(x(1)+h(t,f,h0)-x(4))*(x(1)+h(t,f,h0));
                          g2_best*(x(1)+h(t,f,h0)-x(4))*u(t);
                         -x(2)*x(4)+x(3)*u(t)+theta_m*(x(1)+h(t,f,h0)-x(4))];
     
    [t,x] = ode45(odefun,t,[0,0,0,0]);
    [t,xnoise] = ode45(odefun_noise,t,[0,0,0,0]);
    
    x_pred = xnoise(:,4);
    a_pred = xnoise(:,2);
    b_pred = xnoise(:,3);
    mas = (sum(abs(x(:,1)- x(:,4))))/length(x); %mean absolute error
    mas_noise = (sum(abs(xnoise(:,1)- xnoise(:,4))))/length(x); %mean absolute error when x is recieved with noise
    
    %x and x_pred with and without noise
    fig=figure();
    subplot(2,1,1)
    hold on;
    plot(t,x(:,1));
    plot(t,x(:,4));
    hold off
    grid on;
    xticks(0:5:60)
    title('[Mixed structure] $x$ and $\hat{x}$ without noise','interpreter','latex','FontSize',16);
    xlabel('Time [s]','FontSize',15);
    legend('$x$','$\hat{x}$','interpreter','latex','FontSize',16);
    
    subplot(2,1,2)
    hold on;
    plot(t,xnoise(:,1));
    plot(t,xnoise(:,4));
    hold off
    grid on;
    xticks(0:5:60)
    title({'[Mixed structure] $x$ and $\hat{x}$ with noise';['f =',num2str(f),', h0 =',num2str(h0)]},'interpreter','latex','FontSize',16);
    xlabel('Time [s]','FontSize',15);
    legend('$x$','$\hat{x}$','interpreter','latex','FontSize',16);
    figname="mx f="+f+ "_h0="+ h0+"_x&xhat.jpg";
    saveas(fig,figname)
    
                    
    % e = x - x_pred with and without noise
    fig=figure();
    subplot(2,1,1)
    plot(t,x(:,1)- x(:,4));
    title('[Mixed structure] e = $x$ - $\hat{x}$ without noise','interpreter','latex','FontSize',16);
    xlabel('Time [s]','FontSize',15);
    grid on;
    xticks(0:5:60)
    dim = [.55 .6 .35 .07];
    str = strcat('Mean Absolute Error =',num2str(mas));
    annotation('textbox',dim,'String',str,'FontSize',10)
    
    subplot(2,1,2)
    plot(t,xnoise(:,1)- xnoise(:,4));
    title({'[Mixed structure] e = $x$ - $\hat{x}$ with noise';['f =',num2str(f),', h0 =',num2str(h0)]},'interpreter','latex','FontSize',16);
    xlabel('Time [s]','FontSize',15);
    grid on;
    xticks(0:5:60)
    dim = [.55 .125 .35 .07];
    str = strcat('Mean Absolute Error =',num2str(mas_noise));
    annotation('textbox',dim,'String',str,'FontSize',10)
    figname="mx f="+f+ "_h0="+ h0+"_x-xhat.jpg";
    saveas(fig,figname)
    
    
    
    % a and b with and without noise
    fig=figure();
    subplot(2,1,1)
    hold on;
    plot(t,xnoise(:,2));
    plot(t,x(:,2));
    yline(1.5,'-k');
    hold off
    title({'[Mixed structure] a and $\hat{a}$ with and without noise';['f =',num2str(f),', h0 =',num2str(h0)]},'interpreter','latex','FontSize',16);
    legend('$\hat{a}$ with noise','$\hat{a}$ without noise','$a_{real}$','interpreter','latex');
    xlabel('Time [s]','FontSize',15);
    grid on;
    xticks(0:5:60)
    
    subplot(2,1,2)
    hold on;
    plot(t,xnoise(:,3));
    plot(t,x(:,3));
    yline(2,'-k');
    hold off
    title({'[Mixed structure] b and $\hat{b}$ with and without noise';['f =',num2str(f),', h0 =',num2str(h0)]},'interpreter','latex','FontSize',16);
    legend('$\hat{b}$ with noise','$\hat{b}$ without noise','$b_{real}$','interpreter','latex');
    xlabel('Time [s]','FontSize',15);
    grid on;
    xticks(0:5:60)
    figname="mx f="+f+ "_h0="+ h0+"_a&b.jpg";
    saveas(fig,figname)
    
    fprintf("[Mixed Structure] f = %d and h = %d: Mean absolute error = %f\n",f,h0,mas_noise);
end

