clear, close all
%{
    2022/06/23  bayashi        Lorenz system

    https://en.wikipedia.org/wiki/Lorenz_system
%}
%% Configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 時間積分構成 %%%
Sim.tf        =  100;                          % 終了時刻
Sim.NN        =  100000;                        % 時間積分回数
Sim.dt        =  Sim.tf/Sim.NN;               % 微小時間
Sim.xx0       =  [10; 20; 10];                % 状態変数初期条件
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% システム構成 %%%
Sim.FuncMain  =  @func1;                      % 状態方程式
Sim.val.xx    =  zeros(3, Sim.NN + 1);        % 状態変数
Sim.val.sig   =  10;                          % 
Sim.val.rho   =  28;                          % 
Sim.val.beta  =  8/3;                         % 
Sim.val.A     =  [-Sim.val.sig, Sim.val.sig,     0;
                   Sim.val.rho,          -1,     0;
                     0,  0,  -Sim.val.beta];   % システム行列

%% Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 状態変数時間積分
XX  =  Runge_Kutta(Sim);
xx  =  XX.xx;
% timeaxis
tt = timeMAKE(Sim.NN, Sim.dt);

%% Draw %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for number = 1 : 1
    i = number;
    plot3(xx(1,(1+(i-1)*1E+4):i*1E+4),...
          xx(2,(1+(i-1)*1E+4):i*1E+4),...
          xx(3,(1+(i-1)*1E+4):i*1E+4),...
        'Color',rand(1,3),'LineWidth',0.75);
end
hold on
for number = 2 : 10
    i = number;
    plot3(xx(1,(1+(i-1)*1E+4):i*1E+4),...
          xx(2,(1+(i-1)*1E+4):i*1E+4),...
          xx(3,(1+(i-1)*1E+4):i*1E+4),...
        'Color',rand(1,3),'LineWidth',0.75);
end
hold off
xlabel('X','FontSize',12)
ylabel('Y','FontSize',12)
zlabel('Z','FontSize',12)
grid on
% xlim([-.3 .7])
% ylim([-.4 .3])
% saveas(gcf,'Barchart.png')
exportgraphics(gcf, 'photo_lorenz.png')

%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ss = func1(xx, Sim)
    A  =  Sim.val.A;
    B  =  [0; -xx(1)*xx(3); xx(1)*xx(2)];
    ss =  (A*xx + B).*Sim.dt;
end
function ss = Runge_Kutta(Sim)
    func  =  Sim.FuncMain;
    NN    =  Sim.NN;
    xx    =  Sim.val.xx;
    xx(:,1) =  Sim.xx0;
    for n = 1 : NN
        x0 = xx(:, n);
        k1 = func(x0, Sim);
        k2 = func(x0 + 0.5*k1, Sim);
        k3 = func(x0 + 0.5*k2, Sim);
        k4 = func(x0 + k3, Sim);
        xx(:, n + 1) =  x0 + (k1+2*(k2+k3)+k4)/6;
    end
    ss.xx = xx;
end
function timeaxis  = timeMAKE(Nstep,dt)
    timeaxis = zeros(1,Nstep+1);
    for i = 1:Nstep
        timeaxis(i+1) = dt*i;
    end
end




