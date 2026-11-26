%Filippos Kanlis 9829
%thema 1

close all;
clear;

%System Parameters
m = 10;
b = 0.5;
k = 2.5;
u = @(t) 15*sin(3*t)+8;
tspan= 0:0.1:10;

%solve the Differential equation to get Y values
y0=[0 0];
[t,y] = ode45(@(t,x) fun1(t, x, m, b, k, u), tspan, y0);
y = y(:,1)';
%prepare results of u over tspan
u_t=u(tspan');

%get z values

mhats = [];
bhats = [];
khats = [];
maxPole = 60;

bestError = 100;

mhat = 0;
bhat = 0;
khat = 0;
Bhat = 0;
Khat = 0;
bestPole = 0;
for pole=0:0.5:maxPole
    [mhattemp, bhattemp, khattemp, Bhattemp, Khattemp] = doEstimate(pole, tspan, y, size(tspan,2), u_t);
    mhats(end+1) = mhattemp;
    bhats(end+1) = bhattemp;
    khats(end+1) = khattemp;
    
    error = [m-mhattemp b-bhattemp k-khattemp];
    rmse = sqrt(mean(error.^2));
    
    if rmse < bestError
        mhat = mhattemp;
        bhat = bhattemp;
        khat = khattemp;
        Bhat = Bhattemp;
        Khat = Khattemp;
        bestError = rmse;
        bestPole = pole;
    end
        
end

fprintf("Estimated m: %.3f kg,   Given m: %.3f kg\n", mhat, m);
fprintf("Estimated b: %.3f kg/s,  Given b: %.3f kg/s\n", bhat, b);
fprintf("Estimated k: %.3f kg/s², Given k: %.3f kg/s²\n", khat, k);
fprintf("Best Pole: %f", bestPole);
[t,yhat] = ode45(@(t,x) fun1(t, x, mhat, bhat, khat, u), tspan, y0);
yhat = yhat(:,1)';
error = y - yhat;


%Plots

x = 0:0.5:maxPole;
figYpos=figure(1);
plot(tspan, y)
ylabel("Position y (m)");
xlabel("Time t (s)");
title("y over time");
saveas(figYpos,'figYpos.png');

figEstYpos=figure(2);
plot(tspan, yhat)
ylabel("Estimated Position y-hat");
xlabel("Time t (s)");
title("y-hat over time");
saveas(figEstYpos,'figEstYpos.png');

figEst_m=figure(3);
plot(x, mhats)
yline(m, "Color", 'r', 'LineWidth', 2)
ylabel("m-hat");
xlabel("Pole");
title("m-hat for pole-values");
saveas(figEst_m,'figEst_m.png');

figEst_b=figure(4);
plot(x,bhats)
hold on 
yline(b, "Color", 'r', 'LineWidth', 2)
ylabel("b-hat");
xlabel("Pole");
title("b-hat for pole-values");
saveas(figEst_b,'figEst_b.png');

figEst_k=figure(5);
plot(x,khats)
hold on 
yline(k, "Color", 'r', 'LineWidth', 2)
ylabel("k-hat");
xlabel("Pole");
title("k-hat for pole-values");
saveas(figEst_k,'figEst_k.png');

%Functions
function dx = fun1(t, x, m, b, k, u)

    dx = zeros(2,1);
    dx(1) = x(2);
    dx(2) = 1/m * u(t) - b/m*x(2) - k/m*x(1);
end

function [mhat, bhat, khat, Bhat, Khat] = doEstimate(pole, tspan, y, N, u)
    flt = [1 2*pole pole^2];
    s1 = tf([-1 0],flt);
    s2 = tf(-1,flt);
    s3 = tf(1, flt);

    z1 = lsim(s1,y,tspan);
    z2 = lsim(s2,y,tspan);
    z3 = lsim(s3,u,tspan);

    z = [z1 z2 z3]';

    %Theta calculation
    zSum = zeros(3,1);
    ySum = zeros(3,1);

    for i=1:N
        ztemp = z(:,i);
        ytemp = y(i);

        zSum = zSum + ztemp*ztemp';
        ySum = ySum + ztemp*ytemp;
    end
    theta = ((1/N)*zSum)\((1/N)*ySum);

    %Parameter estimation
    Bhat = theta(1) + 2*pole;
    Khat = theta(2) + pole^2;

    mhat = 1/theta(3);
    bhat = Bhat * mhat;
    khat = Khat * mhat;
end