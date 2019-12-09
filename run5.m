% run 5
% for two coupled neurons, computes the mean interspike interval when
% varying 2 parameters: the coupling strength and the amount of noise added
% to the system
clear all; close all; clc; 

% subunit noise
t = [0:0.01:300];
noise = 'Subunit';
ntrials = 10;
n_neurons = 2;
coupling_strength = 0:0.1:1; 
noise_intensity = 0:0.3:3;
[isi_matrix_su, isi_matrix2_su] = isi_matrix(t, noise, ntrials, n_neurons, noise_intensity, coupling_strength);

% conductance noise
t = [0:0.01:300];
noise = 'FoxLuSystemSize';
ntrials = 10;
n_neurons = 2;
coupling_strength = 0:0.1:1; 
noise_intensity = 0:0.3:3;
[isi_matrix_cd, isi_matrix2_cd] = isi_matrix(t, noise, ntrials, n_neurons, noise_intensity, coupling_strength);

% voltage clamp noise
t = [0:0.01:300];
noise = 'VClamp';
ntrials = 10;
n_neurons = 2;
coupling_strength = 0:0.1:1; 
noise_intensity = 0:0.3:3;
[isi_matrix_vc, isi_matrix2_vc] = isi_matrix(t, noise, ntrials, n_neurons, noise_intensity, coupling_strength);

% markov chain noise
t = [0:0.01:300];
noise = 'MarkovChain';
ntrials = 10;
n_neurons = 2;
coupling_strength = 0:0.1:1; 
noise_intensity = 0:0.3:3;
[isi_matrix_mc, isi_matrix2_mc] = isi_matrix(t, noise, ntrials, n_neurons, noise_intensity, coupling_strength);

save('run5_matrix', 'isi_matrix_su', 'isi_matrix2_su', 'isi_matrix_cd', 'isi_matrix2_cd', ...
    'isi_matrix_vc', 'isi_matrix2_vc', 'isi_matrix_mc', 'isi_matrix2_mc');
