% constructs a matrix of mean interspike intervals
% for varying noise intensities and coupling strengths for three 
% noise models.

% load variables from simulation
load('run5_matrix_subunit', 'isi_matrix_su', 'isi_matrix2_su',
     'isi_matrix_cd', 'isi_matrix2_cd', 'isi_matrix_vc', ...
     'isi_matrix2_vc'); 
load('run5_matrix_new', 'isi_matrix_cd', ...
     'isi_matrix2_cd', 'isi_matrix_vc', 'isi_matrix2_vc');
 
t = [0:0.01:300];
ntrials = 1;
n_neurons = 2;
coupling_strength = 0:0.1:1; 
noise_intensity = 0:0.3:3;

% normalize bounds of colorbar for all of the noise models
bottom = min([min(isi_matrix_su(:)), min(isi_matrix2_su(:)), min(isi_matrix_cd(:)), ... 
    min(isi_matrix2_cd(:)), min(isi_matrix_vc(:)), min(isi_matrix2_vc(:))]);
top = max([max(isi_matrix_su(:)), max(isi_matrix2_su(:)), max(isi_matrix_su(:)),...
    max(isi_matrix2_su(:)), max(isi_matrix_vc(:)), max(isi_matrix2_vc(:))]);
 
% interspike intervals of two coupled neurons for subunit noise 
% first neuron
figure(1);
sgtitle('Mean Interspike Intervals for Subunit Noise', 'fontweight', 'bold', 'fontsize', 15, 'fontweight', 'normal');
subplot(1,2,1);
imagesc(flipud(isi_matrix_su));
caxis([bottom top]);
colorbar;
h = colorbar;
ylabel(h, 'Mean Interspike Interval, <T>', 'fontsize', 15, 'fontweight', 'normal');
title('Neuron 1', 'fontsize', 15, 'fontweight', 'normal');
xlabel('Noise intensity', 'fontsize', 15);
ylabel('Coupling Strength, \kappa', 'fontsize', 15);
c_s = flipud(coupling_strength');
yticks(1:2:length(c_s));
yticklabels(c_s(1:2:end));
xticks(1:2:length(noise_intensity));
xticklabels(noise_intensity(1:2:end));
% second neuron
subplot(1,2,2);
imagesc(flipud(isi_matrix2_su));
caxis([bottom top]);
colorbar; 
h = colorbar;
ylabel(h, 'Mean Interspike Interval, <T>', 'fontsize', 15, 'fontweight', 'normal');
title('Neuron 2','fontsize', 15, 'fontweight', 'normal');
xlabel('Noise intensity','fontsize', 15);
ylabel('Coupling Strength, \kappa','fontsize', 15);
c_s = flipud(coupling_strength');
yticks(1:2:length(c_s));
yticklabels(c_s(1:2:end));
xticks(1:2:length(noise_intensity));
xticklabels(noise_intensity(1:2:end));
set(gcf, 'Position',  [100, 100, 800, 325]);


% interspike intervals of two coupled neurons for system size
% conductance noise
% first neuron
figure(2);
sgtitle('Mean Interspike Intervals for System Size Conductance Noise', 'fontweight', 'bold', 'fontsize', 15, 'fontweight', 'normal');
subplot(1,2,1);
imagesc(flipud(isi_matrix_cd));
caxis([bottom top]);
colorbar;
h = colorbar;
ylabel(h, 'Mean Interspike Interval, <T>', 'fontsize', 15, 'fontweight', 'normal');
title('Neuron 1', 'fontsize', 15, 'fontweight', 'normal');
xlabel('Noise intensity', 'fontsize', 15);
ylabel('Coupling Strength, \kappa', 'fontsize', 15);
c_s = flipud(coupling_strength');
yticks(1:2:length(c_s));
yticklabels(c_s(1:2:end));
xticks(1:2:length(noise_intensity));
xticklabels(noise_intensity(1:2:end));
% second neuron
subplot(1,2,2);
imagesc(flipud(isi_matrix2_cd));
caxis([bottom top]);
colorbar; 
h = colorbar;
ylabel(h, 'Mean Interspike Interval, <T>', 'fontsize', 15, 'fontweight', 'normal');
title('Neuron 2','fontsize', 15, 'fontweight', 'normal');
xlabel('Noise intensity','fontsize', 15);
ylabel('Coupling Strength, \kappa','fontsize', 15);
c_s = flipud(coupling_strength');
yticks(1:2:length(c_s));
yticklabels(c_s(1:2:end));
xticks(1:2:length(noise_intensity));
xticklabels(noise_intensity(1:2:end));
set(gcf, 'Position',  [100, 100, 800, 325]);

% interspike intervals of two coupled neurons for voltage 
% clamp conductance noise
% first neuron
figure(2);
sgtitle('Mean Interspike Intervals for Voltage Clamp Conductance Noise', 'fontweight', 'bold', 'fontsize', 15, 'fontweight', 'normal');
subplot(1,2,1);
imagesc(flipud(isi_matrix_vc));
caxis([bottom top]);
colorbar;
h = colorbar;
ylabel(h, 'Mean Interspike Interval, <T>', 'fontsize', 15, 'fontweight', 'normal');
title('Neuron 1', 'fontsize', 15, 'fontweight', 'normal');
xlabel('Noise intensity', 'fontsize', 15);
ylabel('Coupling Strength, \kappa', 'fontsize', 15);
c_s = flipud(coupling_strength');
yticks(1:2:length(c_s));
yticklabels(c_s(1:2:end));
xticks(1:2:length(noise_intensity));
xticklabels(noise_intensity(1:2:end));
% second neuron
subplot(1,2,2);
imagesc(flipud(isi_matrix2_vc));
caxis([bottom top]);
colorbar; 
h = colorbar;
ylabel(h, 'Mean Interspike Interval, <T>', 'fontsize', 15, 'fontweight', 'normal');
title('Neuron 2','fontsize', 15, 'fontweight', 'normal');
xlabel('Noise intensity','fontsize', 15);
ylabel('Coupling Strength, \kappa','fontsize', 15);
c_s = flipud(coupling_strength');
yticks(1:2:length(c_s));
yticklabels(c_s(1:2:end));
xticks(1:2:length(noise_intensity));
xticklabels(noise_intensity(1:2:end));
set(gcf, 'Position',  [100, 100, 800, 325]);