% plots mean interspike intervals for different noise intensities 
% for subunit, voltage clamp, and system size 
load('run4_isis_witherr', 'isi_vec_su', 'isi_vec_cd', 'isi_vec_vc',... 
    'isi_err_su', 'isi_err_cd', 'isi_err_vc');

figure(1);
hold on;
noise_intensity = 0:0.2:7;
noise_legend = {'Subunit', 'V. Clamp', 'Syst. Size'};
plot(noise_intensity, isi_vec_su, 'b');
plot(noise_intensity, isi_vec_vc, 'color', [0, 0.5, 0.5]);
plot(noise_intensity, isi_vec_cd, 'r')
errorbar(noise_intensity(1:3:end), isi_vec_su(1:3:end), isi_err_su(1:3:end)/sqrt(10), 'b', 'linestyle', 'none');
errorbar(noise_intensity(1:3:end), isi_vec_vc(1:3:end), isi_err_vc(1:3:end)/sqrt(10), 'color', [0 0.5 0.5], 'linestyle', 'none');
errorbar(noise_intensity(1:3:end), isi_vec_cd(1:3:end), isi_err_cd(1:3:end)/sqrt(10), 'r', 'linestyle', 'none');
xlabel('Noise Intensity', 'fontsize', 15); 
ylabel('Mean Interspike Interval, <T>', 'fontsize', 15);
title('Mean Interspike Interval for Varying Noise Intensities', 'fontsize', 15, 'fontweight', 'normal');
legend(noise_legend);
legend boxoff;