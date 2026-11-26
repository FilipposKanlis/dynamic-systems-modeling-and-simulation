% Filippos Kanlis 9829
% Ergasia 2
% thema 1

close all
clc
clear

% erwthma='a' gia u = 10
% erwthma='b' gia u = 10sin(3t) 
erwthma = 'b';

%Meros = 1: Entopismos veltistou g
%Meros = 2: Entopismos veltistou am
%Meros = 3: Methodos Klishs
Meros = 3;

a = 3;
b = 0.5;
fignum = 1;

if erwthma == 'a'
%g=30, am=3 for (a)
    t = 0:0.1:5;
    g = 30;
    am = 3;
    u = @(t) 10;
else
%g=50, am=3 for (b)
    t = 0:0.1:30;
    g = 50;
    am = 3;
    u = @(t) 10*sin(3*t);
end

if Meros == 1
    am=2.8;
    if erwthma =='a'
        g=[1,5,10,20,30,40,50]; %veltisto sto 30
    end
    if erwthma =='b'
        g=[1,5,10,20,30,40,50,60,70]; %veltisto sto 50
    end
    for gi=g 
        odefun = @(t,x) [-a*x(1)+b*u(t);
                        gi*(x(1)-(x(2)*x(4)+x(3)*x(5)))*x(4);
                        gi*(x(1)-(x(2)*x(4)+x(3)*x(5)))*x(5);
                        -am*x(4)+x(1);
                        -am*x(5)+u(t);
                        -am*x(6) + x(2)*x(6) + x(3)*u(t)];
        [t,x] = ode45(odefun,t,[0,0,0,0,0,0]);
        
        x_real = x(:,1);
        x_pred = x(:,6);
        e = x_real - x_pred;
        
        
        fig=figure(fignum);
        fignum = fignum + 1;
        plot(t,e, 'LineWidth', 1, 'Color', 'red');
        xlabel("Time (s)",'Interpreter','Latex','FontSize',16);
        ylabel("Error",'Interpreter','Latex','FontSize',16);
        title(['e = x - $\hat{x}$ for gamma=', num2str(gi)],'interpreter','latex','FontSize',15);
        figname='est_g='+string(gi)+'_part='+erwthma+'.jpg';
        saveas(fig,figname)
    end
end

if Meros == 2
    if erwthma =='a'
        g=30; %veltisto sto 30
        %am=2:0.25:4; %am peripou 3 einai h veltisth timh me g=30
        %am=2.8:0.1:3.3; %am=3 h veltisth timh me g=30    
        am=2.9:0.025:3.1; %am=3 h veltisth timh me g=30
    end
    if erwthma =='b'
        g=50; %veltisto sto 50
        %am=2:0.5:7; %am peripou 3 einai h veltisth timh me g=50
        %am=2.8:0.1:3.3; %am=3 h veltisth timh me g=50    
        am=2.9:0.025:3.1; %am=3 h veltisth timh me g=50
    end

    for ami=am 
        odefun = @(t,x) [-a*x(1)+b*u(t);
                        g*(x(1)-(x(2)*x(4)+x(3)*x(5)))*x(4);
                        g*(x(1)-(x(2)*x(4)+x(3)*x(5)))*x(5);
                        -ami*x(4)+x(1);
                        -ami*x(5)+u(t);
                        -ami*x(6) + x(2)*x(6) + x(3)*u(t)];
        [t,x] = ode45(odefun,t,[0,0,0,0,0,0]);
        
        a_pred = ami-x(:,2);
        b_pred = x(:,3);
        
        
        fig=figure(fignum);
        hold on
        fignum = fignum + 1;
        plot(t,a_pred);
        plot(t,b_pred,'-r');
        yline(a,'--b');
        yline(b,'--r');
        hold off;
        xlabel("Time (s)",'Interpreter','Latex','FontSize',16);
        ylabel("Error",'Interpreter','Latex','FontSize',16);
        title(['$\hat{a}$ and $\hat{b}$ for $a_{m}$=', num2str(ami)],'interpreter','latex','FontSize',16);
        figname='est_am='+string(ami)+'_part='+erwthma+'.jpg';
        saveas(fig,figname)
    end
end

if Meros == 3
    odefun = @(t,x) [-a*x(1)+b*u(t);
                    g*(x(1)-(x(2)*x(4)+x(3)*x(5)))*x(4);
                    g*(x(1)-(x(2)*x(4)+x(3)*x(5)))*x(5);
                    -am*x(4)+x(1);
                    -am*x(5)+u(t);
                    -am*x(6) + x(2)*x(6) + x(3)*u(t)];
    [t,x] = ode45(odefun,t,[0,0,0,0,0,0]);
    
    a_pred = am-x(:,2);
    b_pred = x(:,3);
    x_real = x(:,1);
    x_pred = x(:,6);
    e = x_real - x_pred;

    fig=figure(fignum);
    hold on
    fignum = fignum + 1;
    plot(t,x_pred);
    plot(t,x_real);
    hold off;
    xlabel("Time (s)",'Interpreter','Latex','FontSize',16);
    title(['Plot of $x$ and $\hat{x}$'],'interpreter','latex','FontSize',16);
    figname='x_vs_x_est part='+string(erwthma)+'.jpg';
    saveas(fig,figname)
    
    fig=figure(fignum);
    hold on
    fignum = fignum + 1;
    plot(t,a_pred);
    plot(t,b_pred,'-r');
    yline(a,'--b');
    yline(b,'--r');
    hold off;
    xlabel("Time (s)",'Interpreter','Latex','FontSize',16);
    title('$\hat{a}$ and $\hat{b}$','interpreter','latex','FontSize',16);
    figname='a&b_est_vs_a&b part='+string(erwthma)+'.jpg';
    saveas(fig,figname)

    fig=figure(fignum);
    fignum = fignum + 1;
    plot(t,e, 'LineWidth', 1, 'Color', 'red');
    xlabel("Time (s)",'Interpreter','Latex','FontSize',16);
    ylabel("Error",'Interpreter','Latex','FontSize',16);
    title("Error x-$$\hat{x}$$ ",...
        'Interpreter','Latex','FontSize',16);
    figname='error_x_vs_xest_part='+string(erwthma)+'.jpg';
    saveas(fig,figname)
end