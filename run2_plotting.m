% code for reproducing figure 2 in Goldwyn and Shea-Brown
% plots mean interspike intervals for each noise model for 
% varying levels of input current. includes errorbars to 
% indicate the standard deviation across each of the 10 trials.
load('run2_cv', 'y', 'err', 'cvs', 'cvs_err');

DC_current = 0:1:12;
ntrials = 10;
ntrials_mc = 4;
t = [0:0.01:15000];
noise = {'Markov Chain', 'Subunit', 'VClamp', 'FoxLuSystemSize', 'Current'};
noise_legend = {'Markov Chain', 'Subunit', 'V. Clamp', 'Syst. Size', 'Current'};

% reproducing figure 2A
figure(1);
hold on;
colors = {'k', 'b', [0 0.5 0.5], 'r', 'c'};

plot(DC_current, y(1, :), 'color', 'k');
plot(DC_current, y(2, :), 'color', 'b');
plot(DC_current, y(3, :), 'color', [0 0.5 0.5]);
plot(DC_current, y(4, :), 'color', 'r');
plot(DC_current, y(5, :), 'color', 'c');
plot(DC_current, y(1, :), 'color', 'k');

errorbar(DC_current, y(1, :), err(1, :)/sqrt(4), 'color', 'k');
errorbar(DC_current, y(2, :), err(2, :)/sqrt(10), 'color', 'b');
errorbar(DC_current, y(3, :), err(3, :)/sqrt(10), 'color', [0 0.5 0.5]);
errorbar(DC_current, y(4, :), err(4, :)/sqrt(10), 'color', 'r');
errorbar(DC_current, y(5, :), err(5, :)/sqrt(10), 'color', 'c');
errorbar(DC_current, y(1, :), err(1, :)/sqrt(4), 'color', 'k');

title('Mean Interspike Intervals', 'fontsize', 18, 'fontweight', 'normal');
ylim([0,100])
ylabel('Mean (ms)', 'fontsize', 18);
xlabel('I_{DC} (\muA / cm^2)', 'fontsize', 18);
legend(noise_legend, 'fontsize', 11);
legend boxoff;

%% reproducing figure 2B
figure(2);
hold on;

plot(DC_current, cvs(1, :), 'color', 'k');
plot(DC_current, cvs(2, :), 'color', 'b');
plot(DC_current, cvs(3, :), 'color', [0 0.5 0.5]);
plot(DC_current, cvs(4, :), 'color', 'r');
plot(DC_current, cvs(5, :), 'color', 'c');

errorbar(DC_current, cvs(1, :), cvs_err(1, :)/sqrt(4), 'color', 'k');
errorbar(DC_current, cvs(2, :), cvs_err(2, :)/sqrt(10), 'color', 'b');
errorbar(DC_current, cvs(3, :), cvs_err(3, :)/sqrt(10), 'color', [0 0.5 0.5]);
errorbar(DC_current, cvs(4, :), cvs_err(4, :)/sqrt(10), 'color', 'r');
errorbar(DC_current, cvs(5, :), cvs_err(5, :)/sqrt(10), 'color', 'c');

title('Coefficients of Variance of Interspike Intervals', 'fontsize', 18, 'fontweight', 'normal');
ylabel('CV', 'fontsize', 18);
xlabel('I_{DC} (\mu A / cm^2)', 'fontsize', 18);
legend(noise_legend, 'fontsize', 11);
legend boxoff;
