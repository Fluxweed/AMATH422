%testbench for noised coupled neurons
close all
clear all
clc

% without coupling terms 
k = 0.0
varSim = noised_coupling(5,k, [0:0.01:200], @(t)sin(t), [], 10, 'Subunit');
t = varSim.t';
cmap = colormap(lines);

for i =1:size(varSim.V,2)
    plot(t, varSim.V(:,i),'Color', cmap(i,:)); hold on;
end
%
k = 0.9
%varSim = noised_coupling(5,k, [0:0.01:200], @(t)sin(t), [], 10, 'FoxLuSystemSize');
t = varSim.t';
cmap = colormap(lines);
figure
for i =1:size(varSim.V,2)
    plot(t, varSim.V(:,i),'Color', cmap(i,:)); hold on;
end




%get variables

n1= varSim.V(:,1);
n2 =varSim.V(:,2);
%n3 = varSim.V(:,3);



% plot(t, n3, 'b');
%figure(2)
%plot(t, varSim.I)