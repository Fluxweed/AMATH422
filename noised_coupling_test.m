%testbench for noised coupled neurons
close all
clear all

for k = 0:0.1:0.5
    [Y varSim] = noised_coupling(2,5,k, [0:0.01:200], @(t)sin(t), [], 10, 'Subunit');
    t = varSim.t';
    cmap = colormap(lines);
    figure
    for i =1:size(varSim.V,2)
        plot(t, varSim.V(:,i),'Color', cmap(i,:)); hold on;
    end
end


%get variables

n1= varSim.V(:,1);
n2 =varSim.V(:,2);
%n3 = varSim.V(:,3);



% plot(t, n3, 'b');
%figure(2)
%plot(t, varSim.I)