%Filippos Kanlis 9829
%thema 2

close all;
clear;


% VC and VR values
tspan = 0:10^-4:100;
[vr, vc] = v(tspan);
vr = vr.';
vc = vc.';


N = length(tspan);

pole = 30;


% Issue A
  % 1 to not search for the optimal pole. 1 to search.
    % Defaulted at 1.
optimize = 1;
    
if optimize
    [invRC, invLC, pole] = findOptimalEstimates(vc, N, tspan, 200);
    fprintf("Best pole at: %f\n", pole);
else
    [invRC, invLC, mse] = findEstimates(vc, N, tspan, pole);
end
makeErrorPlots(vc, vr, invRC, invLC, tspan);

qaInvRC = 0.1503;
qaInvLC = 1000600;

% Issue B
randElements = 5;
randIdx = randperm(N,randElements);

for r=1:randElements
    ridx = randIdx(r);
    vc(ridx) = vc(ridx) + 10^2;
end
maxPole = 400;

rcerr = [];
lcerr = [];

for p = 0:10:maxPole
    [invRC, invLC, mser] = findEstimates(vc, N, tspan, p);
    rcerr(end+1) = qaInvRC - invRC;
    lcerr(end+1) = qaInvLC - invLC;
end
x = [0:10:maxPole];
fig5=figure(5);
plot(x, rcerr, 'Color', '#A2142F')  
ylabel("$$\frac{1}{RC}$$ error",'Interpreter','Latex');
xlabel("Pole",'Interpreter','Latex');
title("$$\hat{\frac{1}{RC}} - \hat{\frac{1}{RC}}_B$$ for different poles",'Interpreter','Latex');
saveas(fig5,'fig5.png');

fig6=figure(6);
plot(x,lcerr)
ylabel("$$\frac{1}{LC}$$ error",'Interpreter','Latex');
xlabel("Pole");
title("$$\hat{\frac{1}{LC}} - \hat{\frac{1}{LC}}_B$$ for different poles",'Interpreter','Latex');
saveas(fig6,'fig6.png');

%findOptimalEstimates(vc, N, tspan, maxPole);
    


% Functions
function [invRC, invLC, mser] = findEstimates(vc, N, tspan, pole)

    % Get z values
    filter = [1 2*pole pole^2];
    sys1 = tf([-1 0],filter);
    sys2 = tf([-1],filter);
    sys4 = tf([1 0], filter);
    sys3 = tf([1], filter);

    u1 = 2*sin(4*tspan);
    u2 = zeros(1, N) + 4;

    z1 = lsim(sys1,vc,tspan);
    z2 = lsim(sys2,vc,tspan);
    z3 = lsim(sys3,u1,tspan);
    z4 = lsim(sys4,u1,tspan);
    z5 = lsim(sys3,u2,tspan);
    z6 = lsim(sys4,u2,tspan);

    z = [z1 z2 z3 z4 z5 z6]';

    % Find theta
    zum = zeros(6,1);
    yum = zeros(6,1);

    for i=1:N
        ztemp = z(:,i);
        ytemp = vc(i);

        zum = zum + ztemp*(ztemp');
        yum = yum + ztemp*ytemp;
    end

    theta = ((1/N)*zum)\((1/N)*yum);

    % Retrieve parameters
    theta_star = cat(2,[theta(1) + 2*pole, theta(2) + pole^2], theta(3:end).')';

    % Plot error

    invRC = theta_star(1);
    invLC = theta_star(2);

    vchat = solveVCode(tspan,invLC, invRC);
    %yhat = ((theta')*z)';
    error = vchat - vc;
    %plot(tspan, error)
    mser = mse(vchat, vc);%mean(abs(error));
end

function [invRC, invLC, pole] = findOptimalEstimates(vc, N, tspan, maxPole)
    bestError = 10^2;
    invRC = 0;
    invLC = 0;
    pole = 0;
    invRCArray = [];
    invLCArray = [];
    for p=0:10:maxPole
        [tempinvRC, tempinvLC, mse] = findEstimates(vc, N, tspan, p);
        invRCArray(end+1) = tempinvRC;
        invLCArray(end+1) = tempinvLC;
        if mse < bestError
            bestError = mse;
            invRC = tempinvRC;
            invLC = tempinvLC;
            pole = p;
        end
    end

    fig1=figure(1);
    x = [0:10:maxPole];
    plot(x, invRCArray,"Color","#7E2F8E")
    ylabel("$$\hat{\frac{1}{RC}}$$",'Interpreter','Latex');
    xlabel("Pole");
    title("$$\frac{1}{RC}$$ estimate for different poles",'Interpreter','Latex');
    saveas(fig1,'fig1.png');

    fig2=figure(2);
    plot(x, invLCArray, "Color", '#D95319')
    ylabel("$$\hat{\frac{1}{LC}}$$",'Interpreter','Latex');
    xlabel("Pole");
    title("$$\frac{1}{LC}$$ estimate for different poles",'Interpreter','Latex');
    saveas(fig2,'fig2.png');
end
function yp = odefun(t, y, invLC, invRC)

    yp = zeros(2,1);
    yp(1) = y(2);
    
    u1 = 3*sin(2*t);
    u1dot = 6*cos(2*t);
    
    u2 = 2;
    u2dot = 0;

    yp(2) = -invRC*y(2) - invLC*y(1) + invLC*u2 + invRC*u1dot + invRC*u2dot;
end

function [vc] = solveVCode(tspan,invLC, invRC)
    % initial values
    initStates = [0,0];
    
    % Get VC values
    [t,vc] = ode45(@(t,y) odefun(t, y, invLC, invRC), tspan, initStates);

    vc = vc(:,1);
end

function makeErrorPlots(vc, vr, invRC, invLC, tspan)
    N = length(tspan);
    fprintf("Estimate for 1/RC = %.3f\n",invRC)
    fprintf("Estimate for 1/LC = %.2f",invLC)

    vchat = solveVCode(tspan,invLC, invRC);
    vc_error = vc - vchat;

    u1 = 3*sin(2*tspan);
    u2 = zeros(1, N) + 2;
    vrhat = (u1 + u2)' - vchat;
    vr_error = vr - vrhat;

    fig31=figure(31);
    plot(tspan, vc, "Color", 'r')
    ylabel("Error");
    xlabel("Time (s)");
    title("V_c vs time",'Interpreter','Latex');
    saveas(fig31,'fig31.png');

    fig32=figure(32);
    plot(tspan, vchat, "Color", 'r')
    ylabel("Error");
    xlabel("Time (s)");
    title("$$\hat{V_c}$$ vs time",'Interpreter','Latex');
    saveas(fig32,'fig32.png');

    fig33=figure(33);
    plot(tspan, vc_error, "Color", 'r')
    ylabel("Error");
    xlabel("Time (s)");
    title("$$V_c - \hat{V_c}$$ vs time",'Interpreter','Latex');
    saveas(fig33,'fig33.png');

    fig41=figure(41);
    plot(tspan, vr, "Color", 'g')
    ylabel("Error");
    xlabel("Time (s)");
    title("$$V_r$$ over time",'Interpreter','Latex');
    saveas(fig41,'fig41.png');

    fig42=figure(42);
    plot(tspan, vrhat, "Color", 'g')
    ylabel("Error");
    xlabel("Time (s)");
    title("$$\hat{V_r}$$ over time",'Interpreter','Latex');
    saveas(fig42,'fig42.png');

    fig43=figure(43);
    plot(tspan, vr_error, "Color", 'g')
    ylabel("Error");
    xlabel("Time (s)");
    title("$$V_r - \hat{V_r}$$ over time",'Interpreter','Latex');
    saveas(fig43,'fig43.png');
end