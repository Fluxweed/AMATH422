% for two coupled neurons, computes the mean interspike interval when
% varying 2 parameters: the coupling strength and the amount of noise added
% to the system. plots the mean interspike intervals for different 
% parameter combinations on a colormap.
clear all; close all; clc; 

t = [0:0.01:100];
noise = 'Subunit';
ntrials = 10;
coupling_strength = 0:0.1:1;
noise_intensity = 1:0.2:3; 
isi_matrix = zeros(length(coupling_strength), length(noise_intensity));
isi_matrix2 = zeros(length(coupling_strength), length(noise_intensity));
for i = 1:length(coupling_strength)
    for j = 1:length(noise_intensity)
        kappa = coupling_strength(i);
        n_s = noise_intensity(j);
        [avg_mean_isi, avg_mean_isi2] = mean_isi_coupled(t, noise, ntrials, kappa, n_s);
        isi_matrix(i, j) = avg_mean_isi;
        isi_matrix2(i, j) = avg_mean_isi2;
    end 
end 

bottom = min(min(isi_matrix(:)), min(isi_matrix2(:)));
top = max(max(isi_matrix(:)), max(isi_matrix2(:)));

figure(1);
sgtitle('Mean Interspike Intervals for Subunit Noise');

subplot(1,2,1);
imagesc(flipud(isi_matrix));
caxis([bottom top]);
colorbar;
h = colorbar;
ylabel(h, 'Mean Interspike Interval, <T>', 'fontsize', 11);
title('Neuron 1', 'fontsize', 11);
xlabel('Noise intensity', 'fontsize', 11);
ylabel('Coupling Strength, \kappa', 'fontsize', 11);
c_s = flipud(coupling_strength');
yticks(1:2:length(c_s));
yticklabels(c_s(1:2:end));
xticks(1:2:length(noise_intensity));
xticklabels(noise_intensity(1:2:end));

subplot(1,2,2);
imagesc(flipud(isi_matrix2));
caxis([bottom top]);
colorbar; 
h = colorbar;
ylabel(h, 'Mean Interspike Interval, <T>', 'fontsize', 11);
title('Neuron 2','fontsize', 11);
xlabel('Noise intensity','fontsize', 11);
ylabel('Coupling Strength, \kappa','fontsize', 11);
c_s = flipud(coupling_strength');
yticks(1:2:length(c_s));
yticklabels(c_s(1:2:end));
xticks(1:2:length(noise_intensity));
xticklabels(noise_intensity(1:2:end));

set(gcf, 'Position',  [100, 100, 800, 350]);

%% Difference in Mean Interspike Intervals
% plot the residuals of the mean interspike intervals for neuron 1 and
% neuron 2

figure(2);
dif_matrix = flipud(isi_matrix - isi_matrix2);
imagesc(dif_matrix);
lim = max(abs(min(dif_matrix(:))), abs(max(dif_matrix(:))));
caxis([-lim, lim]);
colorbar;
h = colorbar;
colormap(redblue);
ylabel(h, '<T>_1 - <T>_2', 'fontsize', 11);
title('Difference in Mean Interspike Intervals for Coupled Neurons', 'fontsize', 11);
xlabel('Noise intensity', 'fontsize', 11);
ylabel('Coupling Strength, \kappa', 'fontsize', 11);
c_s = flipud(coupling_strength');
yticks(1:2:length(c_s));
yticklabels(c_s(1:2:end));
xticks(1:2:length(noise_intensity));
xticklabels(noise_intensity(1:2:end));


