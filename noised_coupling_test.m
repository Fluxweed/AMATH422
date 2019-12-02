%testbench for noised coupled neurons
close all
clear all

[Y varSim] = noised_coupling(2,10,0.3, [0:0.01:100], @(t)sin(t), [], 100, 'Subunit');


%get variables
t = varSim.t';
n1= varSim.V(:,1);
n2 =varSim.V(:,2);
%n3 = varSim.V(:,3);


cmap = colormap(lines);
for i =1:size(varSim.V,2)
    plot(t, varSim.V(:,i),'Color', cmap(i,:)); hold on;
end
% plot(t, n3, 'b');
figure(2)
plot(t, varSim.I)