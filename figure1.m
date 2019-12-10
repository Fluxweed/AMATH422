% code to get figure 1
clear all;
close all;
% simulate variance 
nSim = 5000;
sim =1; 
NoiseModel = 'FoxLuSystemSize';
t = [0:0.01:100];
SigmaIn = [];
Area = 100;
NNa = 6000;
NK =1800;

[Y varSim] = simulate_variance(sim, nSim, t, @(t)0,SigmaIn, Area, NoiseModel);


%calculate m h n
m  = Y(:,5);
h = Y(:,6);
n = Y(:,7);


na = varSim.NaFraction;
k = varSim.KFraction;


meanNa_mc =( m.^3).*h;
meanK_mc = n.^4;
meanNa = mean(na,2);
meanK = mean(k,2);
figure(1);
%figure 1.a
subplot(4,1,1);
plot(Y(:,1), Y(:,2));xlabel('Time (ms)');
ylabel('Voltage (mV)');
subplot(4,1,2);
%figure 1.b
plot(t, meanNa_mc, '--'); hold on;
plot(t, meanNa, 'c');
plot(t, meanK, 'g');
plot(t, meanK_mc, '--'); xlabel('Time (ms)');
legend('Na','K');
ylabel('Mean Open Fraction');

%calculate variance
varNa = var(na, [],2);
varK = var(k, [],2);
varNa_mc= (m.^3.*h.*(1-m.^3))./NNa
varK_mc = (n.^4.*(1-n.^4))./NK
subplot(4,1,3);

%figure 1.c
plot(Y(:,1), zscore(varNa), 'r*'); hold on;
plot(Y(:,1), zscore(varNa_mc), 'k','LineWidth', 4); xlabel('Time (ms)');
legend('Syst. Size');
ylabel('Var Open Na Fraction');
subplot(4,1,4);


%figure 1.d
plot(Y(:,1), zscore(varK) ,'r*'); hold on;
plot(Y(:,1), zscore(varK_mc), 'k','LineWidth', 4);  xlabel('Time (ms)');
legend('Syst. Size');
ylabel('Var Open K Fraction');
saveas(gcf, 'f1','epsc');

% normNa = 1000;
% normK = 100000;
% figure('Position',[200 200 1000 900]);
% 
% subplot(4,1,1);
% plot(Y(:,1), Y(:,2));xlabel('Time (ms)');
% ylabel('Voltage (mV)');
% subplot(4,1,2);
% plot(t, meanNa_mc, '--'); hold on;
% plot(t, meanNa, 'c');
% plot(t, meanK, 'g');
% plot(t, meanK_mc, '--'); xlabel('Time (ms)');legend('Na','K');
% 
% ylabel('Mean Open Fraction');
% %calculate variance
% varNa = var(na, [],2);
% varK = var(k, [],2);
% varNa_mc= (m.^3.*h.*(1-m.^3))./NNa
% varK_mc = (n.^4.*(1-n.^4))./NK
% subplot(4,1,3);
% plot(Y(:,1), normNa*varNa, 'r*'); hold on;
% plot(Y(:,1),varNa_mc, 'k','LineWidth', 4); xlabel('Time (ms)');legend('Syst. Size');
% 
% ylabel('Var Open Na Fraction');
% subplot(4,1,4);
% plot(Y(:,1), normK*varK ,'b*'); hold on;
% plot(Y(:,1), varK_mc, 'k','LineWidth', 4);  xlabel('Time (ms)');legend('Syst. Size');
% 
% 
% ylabel('Var Open K Fraction');

saveas(gcf,'figure1', 'epsc');




