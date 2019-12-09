% run2, reproducing figure 2
% compute mean interspike intervals for each noise model, then plot
% for varying amounts of input current
clear all; close all; clc; 

DC_current = 0:1:12;
ntrials = 10;
ntrials_mc = 4;
t = [0:0.01:25000];
noise = string({'Markov Chain', 'Subunit', 'VClamp', 'FoxLuSystemSize', 'Current'});

[mc_mean_isi, mc_std_isi, mc_cv, mc_cv_std] = isi(t, noise(1), DC_current, ntrials_mc); 
[s_mean_isi, s_std_isi, s_cv, s_cv_std] = isi(t, noise(2), DC_current, ntrials); 
[vc_mean_isi, vc_std_isi, vc_cv, vc_cv_std] = isi(t, noise(3), DC_current, ntrials);
[ss_mean_isi, ss_std_isi, ss_cv, ss_cv_std] = isi(t, noise(4), DC_current, ntrials);
[c_mean_isi, c_std_isi, c_cv, c_cv_std] = isi(t, noise(5), DC_current, ntrials); 

y = [mc_mean_isi; s_mean_isi; vc_mean_isi; ss_mean_isi; c_mean_isi];
err = [mc_std_isi; s_std_isi; vc_std_isi; ss_std_isi; c_std_isi];

cvs = [mc_cv; s_cv; vc_cv; ss_cv; c_cv];
cvs_err = [mc_cv_std; s_cv_std; vc_cv_std; ss_cv_std; c_cv_std];


save('run2_cv', 'y', 'err', 'cvs', 'cvs_err');