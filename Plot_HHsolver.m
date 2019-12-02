close all; clear all; clc;
trials = 100;
time = 0:0.01:200;
FractionN = zeros(trials,length(time));
FractionK = zeros(trials,length(time));
VOut = zeros(trials,length(time));
VOut2 = zeros(trials,length(time));
Subunit = zeros(1,trials);
VClamp = zeros(1,trials);
System = zeros(1,trials);

for i = 0:0.1:1;
    %E = StochasticHH_func(time, @(t) 0, [], 100,'ODE');
    
    %S = StochasticHH_func([0:0.01:100], @(t) 0, [1], 100,'Subunit');
    %V = StochasticHH_func([0:0.01:100], @(t) 0, [1], 100,'VClamp');
    FLS = StochasticHH_funcb(time, @(t) 0, [1], 100,'Subunit',i);
    %FractionN(i,:) = FLS(:,3);
    %FractionK(i,:) = FLS(:,4);
    %VOut(i,:) = FLS(:,2);
    %VOut2(i,:) = FLS(:,8);
    
    figure
    hold on;
    plot(FLS(:,1), FLS(:,2))
    plot(FLS(:,1), FLS(:,9))
    title('Voltage vs. Time')
end

hold on;
%figure(1)
%plot(FLS(:,1), FLS(:,3)), title('Fraction of Open Na channels vs. Time')
%figure(2)
%plot(FLS(:,1), FLS(:,4)), title('Fraction of Open K channels vs. Time')
%figure(1)
%plot(FLS(:,1), FLS(:,2)), title('Voltage vs. Time')
%plot(FLS(:,1), FLS(:,9)), title('Voltage vs. Time 2')
%figure(2)
%plot(FLS(:,1), FLS(:,3)), title('Fraction of Open Na channels vs. Time')
%plot(FLS(:,1), FLS(:,10)), title('Fraction of Open Na channels vs. Time 2')
%figure(3)
%plot(FLS(:,1), FLS(:,16)), title('Voltage vs. Time 3')
