% Model parameters
%beta = 0.0460; % rate of infection 0.0366
%gamma = 0.0093; % rate of recovery (try also 0.07)

N = 7545198; % Total population N = S + I + R
I0 = 221.5238; % initial number of infected
T = 365; % period of 300 days
dt = 1; % time interval of 6 hours (1/4 of a day)
%D0 = 302; % death count at day 0
%f = 0.0152; % fatality rate 
R0 = 2165;

%beta = 0.0408;
%gamma = 0.0153;
%f = 0.0166;
e1 = 0.5;
e2 = 0.9;
p1 = 0.5;
p2 = 0.5;
% beta = 0.0379 * (1 - (p1*israelRate(1)*e1 + p2*israelRate(1)*e2));
% beta = table2array(israelrate2);
beta = 0.139708;
% beta = 0.17109;
%beta = 0.2;
beta = 0.184449;
%beta = 0.3615;

gamma = 0.040345;
f = 0.0166;

delta = 1/280; % rate of immunity loss
% Calculate the model
[S,I,R,IC,dIT,R0,D,vaccineRate] = sir_model(beta,gamma,delta,N,I0,T,dt,R0);
% Curve
tt = 0:dt:T-dt;
subplot(1,4,1);
plot(tt,D,'black',tt,I,'r',tt,R,'g',tt,IC,'magenta','LineWidth',2); grid on;
xlim([0,180])
%plot(tt,S(1));
xlabel('Days');ylabel('Number of individuals');
legend('D(death)','I(infected)','R(recovered)','Cummulative Infection');
subplot(1,4,2);
plot(tt,dIT,'r','LineWidth',2); grid on;
xlabel('Days');ylabel('Number of individuals');
legend('New infections per day');
xlim([0,80]);
subplot(1,4,3);
plot(tt,D,'black','LineWidth',2); grid on;
xlim([0,180])
%plot(tt,S(1));
xlabel('Days');ylabel('Number of individuals');
legend('D(death)');
subplot(1,4,4);
plot(tt,vaccineRate,'r',tt,R0,'b', 'LineWidth', 2); grid on;
xlim([0,80])
h = yline(1, 'r--', 'LineWidth', 4);
legend('vaccination %', 'R0');
%plot(R0,'black','LineWidth',2); grid on;
%xlabel('Days');ylabel('R0');
%legend('R0');
%ylim([0,3])
%fprintf('Value of parameter R0 is %.2f',beta/gamma)

function [S,I,R,IC,dIT,R0,D,vaccineRate] = sir_model(beta,gamma,delta,N,I0,T,dt, R0)
    % if delta = 0 we assume a model without immunity loss
    S = zeros(1,T/dt);
    S(1) = N;
    I = zeros(1,T/dt);
    I(1) = I0;
    R = zeros(1,T/dt);
    R(1) = R0;
    D = zeros(1,T/dt);
    D(1) = 0;
    IC = zeros(1, T/dt);
    IC(1) = I0;
    dIT = zeros(1,T/dt);
    temp = 0.00197820202;
    R0 = zeros(1,T/dt);
    f = 0.017870885;
    f
    vaccineRate = zeros(1,T/dt);
    
    for tt = 1:(T/dt)-1
        % Equations of the model
        %if tt <= 218
        %    beta = B(tt);
        %else
        %    beta = B(218);
        %end
        %beta = beta * (1 - temp *tt * 0.6594 * 0.538 - temp *tt* (1-0.6594) * 0.9);
        beta = beta * (1 - temp*tt*0.91);
        D(tt+1) = f * IC(tt);
        dD = D(tt+1) - D(tt);
        
        %dS = (-beta*I(tt)*S(tt)/N + delta*R(tt) -  S(tt)*temp*(tt)*(coronavacRatio * 0.5 + (1 - coronavacRatio) * 0.9 ) - dD);
        dS = (-beta*I(tt)*S(tt)/N + delta*R(tt));
        % dI = (beta*I(tt)*S(tt)/N - gamma*I(tt) - (1-f)*I(tt));
        % dI = ((beta*I(tt)*S(tt)/N - gamma*I(tt));
        % dI = (beta*I(tt)*S(tt)/N - gamma*I(tt));
        dI = (beta*I(tt)*S(tt)/N - gamma*I(tt));
        dR = (gamma*I(tt) - delta*R(tt));
        S(tt+1) = S(tt) + dS;
        I(tt+1) = I(tt) + dI;
        R(tt+1) = R(tt) + dR;
        D(tt+1) = f * IC(tt);
        R0(tt+1) = (-dS*N/I(tt)/S(tt)) / (dR/I(tt));
        IC(tt+1) = IC(tt) + (beta*I(tt)*S(tt)/N);
        dIT(tt) = (beta*I(tt)*S(tt)/N);      
        R0(tt) = beta/(gamma+f);
        vaccineRate(tt) = temp*tt*100;
      
    end
end