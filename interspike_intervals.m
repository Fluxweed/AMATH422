%% reproducing figure 2
% compute mean interspike intervals for each noise model, then plot
% for varying amounts of input current
clear all; close all; clc; 

DC_current = 0:1:12;
ntrials = 10;
ntrials_mc = 4;
t = [0:0.01:1000];
noise = string({'Markov Chain', 'Subunit', 'VClamp', 'FoxLuSystemSize', 'Current'});
noise_legend = string({'Markov Chain', 'Subunit', 'V. Clamp', 'Syst. Size', 'Current'});

[mc_mean_isi, mc_std_isi, mc_cv, mc_cv_std] = isi(t, noise(1), DC_current, ntrials_mc); 
disp('1');
[s_mean_isi, s_std_isi, s_cv, s_cv_std] = isi(t, noise(2), DC_current, ntrials); 
disp('2');
[vc_mean_isi, vc_std_isi, vc_cv, vc_cv_std] = isi(t, noise(3), DC_current, ntrials);
disp('3');
[ss_mean_isi, ss_std_isi, ss_cv, ss_cv_std] = isi(t, noise(4), DC_current, ntrials);
disp('4');
[c_mean_isi, c_std_isi, c_cv, c_cv_std] = isi(t, noise(5), DC_current, ntrials); 

% figure 2A
figure(1);
hold on;
colors = ['k', 'b', 'g', 'r', 'c'];
y = [mc_mean_isi; s_mean_isi; vc_mean_isi; ss_mean_isi; c_mean_isi];
err = [mc_std_isi; s_std_isi; vc_std_isi; ss_std_isi; c_std_isi];

for i = 1:size(y, 1)
    plot(DC_current, y(i, :), colors(i));
end

for i = 1:size(y, 1)
    errorbar(DC_current, y(i, :), err(i, :), colors(i));
end

title('Mean Interspike Intervals');
ylabel('Mean (ms)');
xlabel('I_{DC} (\muA / cm^2)');
legend(noise_legend);

%% figure 2B
figure(2);
hold on;

cvs = [mc_cv; s_cv; vc_cv; ss_cv; c_cv];
cvs_err = [mc_cv_std; s_cv_std; vc_cv_std; ss_cv_std; c_cv_std];

for i = 1:size(cvs, 1)
    plot(DC_current, cvs(i, :), colors(i));
end

for i = 1:size(cvs, 1)
    errorbar(DC_current, cvs(i, :), cvs_err(i, :), colors(i));
end

title('Coefficients of Variance of Interspike Intervals');
ylabel('CV');
xlabel('I_{DC} (\mu A / cm^2)');
legend(noise_legend);
