% for two coupled neurons, computes the mean interspike interval when
% varying 2 parameters: the coupling strength and the amount of noise added
% to the system. plots the mean interspike intervals for different 
% parameter combinations on a colormap.
clear all; close all; clc; 

%% subunit noise
t = [0:0.01:500];
noise = 'Subunit';
ntrials = 1;
n_neurons = 2;
coupling_strength = 0:0.1:1; 
noise_intensity = 0:0.3:3;
[isi_matrix_subunit, isi_matrix2_subunit] = isi_matrix(t, noise, ntrials, n_neurons, noise_intensity, coupling_strength);

bottom = min(min(isi_matrix_subunit(:)), min(isi_matrix2_subunit(:)));
top = max(max(isi_matrix_subunit(:)), max(isi_matrix2_subunit(:)));

figure(1);
sgtitle('Mean Interspike Intervals for Subunit Noise', 'fontweight', 'bold', 'fontsize', 14);

subplot(1,2,1);
imagesc(flipud(isi_matrix_subunit));
caxis([bottom top]);
colorbar;
h = colorbar;
ylabel(h, 'Mean Interspike Interval, <T>', 'fontsize', 13);
title('Neuron 1', 'fontsize', 13);
xlabel('Noise intensity', 'fontsize', 13);
ylabel('Coupling Strength, \kappa', 'fontsize', 13);
c_s = flipud(coupling_strength');
yticks(1:2:length(c_s));
yticklabels(c_s(1:2:end));
xticks(1:2:length(noise_intensity));
xticklabels(noise_intensity(1:2:end));

subplot(1,2,2);
imagesc(flipud(isi_matrix2_subunit));
caxis([bottom top]);
colorbar; 
h = colorbar;
ylabel(h, 'Mean Interspike Interval, <T>', 'fontsize', 13);
title('Neuron 2','fontsize', 13);
xlabel('Noise intensity','fontsize', 13);
ylabel('Coupling Strength, \kappa','fontsize', 13);
c_s = flipud(coupling_strength');
yticks(1:2:length(c_s));
yticklabels(c_s(1:2:end));
xticks(1:2:length(noise_intensity));
xticklabels(noise_intensity(1:2:end));

set(gcf, 'Position',  [100, 100, 800, 325]);

%% conductance noise 
t = [0:0.01:500];
noise = 'FoxLuSystemSize';
ntrials = 1;
n_neurons = 2;
coupling_strength = 0:0.1:1; 
noise_intensity = 0:0.3:3;
isi_matrix_cond = zeros(length(coupling_strength), length(noise_intensity));
isi_matrix2_cond = zeros(length(coupling_strength), length(noise_intensity));
for i = 1:length(coupling_strength)
    for j = 1:length(noise_intensity)
        kappa = coupling_strength(i);
        n_s = noise_intensity(j);
        [avg_mean_isi, avg_mean_isi2] = mean_isi_coupled(t, noise, ntrials, n_neurons, kappa, n_s);
        isi_matrix_cond(i, j) = avg_mean_isi;
        isi_matrix2_cond(i, j) = avg_mean_isi2;
    end 
end 

bottom = min(min(isi_matrix_cond(:)), min(isi_matrix2_cond(:)));
top = max(max(isi_matrix_cond(:)), max(isi_matrix2_cond(:)));

figure(1);
sgtitle('Mean Interspike Intervals for Conductance Noise', 'fontweight', 'bold', 'fontsize', 14);

subplot(1,2,1);
imagesc(flipud(isi_matrix_cond));
caxis([bottom top]);
colorbar;
h = colorbar;
ylabel(h, 'Mean Interspike Interval, <T>', 'fontsize', 13);
title('Neuron 1', 'fontsize', 13);
xlabel('Noise intensity', 'fontsize', 13);
ylabel('Coupling Strength, \kappa', 'fontsize', 13);
c_s = flipud(coupling_strength');
yticks(1:2:length(c_s));
yticklabels(c_s(1:2:end));
xticks(1:2:length(noise_intensity));
xticklabels(noise_intensity(1:2:end));

subplot(1,2,2);
imagesc(flipud(isi_matrix2_cond));
caxis([bottom top]);
colorbar; 
h = colorbar;
ylabel(h, 'Mean Interspike Interval, <T>', 'fontsize', 13);
title('Neuron 2','fontsize', 13);
xlabel('Noise intensity','fontsize', 13);
ylabel('Coupling Strength, \kappa','fontsize', 13);
c_s = flipud(coupling_strength');
yticks(1:2:length(c_s));
yticklabels(c_s(1:2:end));
xticks(1:2:length(noise_intensity));
xticklabels(noise_intensity(1:2:end));

set(gcf, 'Position',  [100, 100, 800, 325]);

%% Difference in Mean Interspike Intervals
% plot the differences between the mean interspike intervals for 
% neuron 1 and neuron 2

figure(2);
dif_matrix = flipud(isi_matrix2_subunit - isi_matrix2_cond);
imagesc(dif_matrix);
lim = max(abs(min(dif_matrix(:))), abs(max(dif_matrix(:))));
caxis([-lim, lim]);
colorbar;
h = colorbar;
colormap(redblue);
ylabel(h, '<T>(subunit) - <T>(conductance)', 'fontsize', 13);
title('Difference in Mean Interspike Intervals for Coupled Neurons', 'fontsize', 13);
xlabel('Noise intensity', 'fontsize', 13);
ylabel('Coupling Strength, \kappa', 'fontsize', 13);
c_s = flipud(coupling_strength');
yticks(1:2:length(c_s));
yticklabels(c_s(1:2:end));
xticks(1:2:length(noise_intensity));
xticklabels(noise_intensity(1:2:end));



%%
figure(3);
hold on;
plot(noise_intensity, isi_vec_subunit, 'k');
plot(noise_intensity, isi_vec_cond, 'r')
xlabel('Noise Intensity', 'fontsize', 16); 
ylabel('Mean Interspike Interval, <T>', 'fontsize', 16);
title('Mean Interspike Interval for Different Types of Noise', 'fontsize', 16);
legend('Subunit', 'Conductance');





