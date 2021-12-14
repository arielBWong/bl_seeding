function test_kriging
clear all;close all;clc;
addpath E:\BING\OODACE
startup;

% Demonstration of cheap and expensive estimates: Cokriging (cheap first)
x_trg{1,:}=[0,0.2,0.4,0.6,0.8,1]';
y_trg{1,:}=eval_fun(x_trg{1});

x_trg{2,:}=[0,0.5,1]';
y_trg{2,:}=eval_fun(x_trg{2});

% Fitting the model
k = oodacefit(x_trg,y_trg);

% Testing the performance
x_tst=[linspace(0,1,100)]';
y_tst=eval_fun(x_tst);
[y_prd,~] = k.predict(x_tst);
plot(x_tst,y_tst,'r-');hold on;
plot(x_tst,y_prd,'b-');
plot(x_trg{1,:},y_trg{1,:},'ro','Markersize',15); hold on;
if(size(x_trg,1)>1)
    plot(x_trg{2,:},y_trg{2,:},'ro','Markerfacecolor','r','Markersize',15); hold on;
    legend('True function','Approximated function','Cheap Samples','Expensive Samples','FontSize',16);
else
    
    plot(x_trg{1,:},y_trg{1,:},'ro','Markersize',15); hold on;
    legend('True function','Approximated function','Cheap Samples','FontSize',16);
end
xlabel('Variable','FontSize',16);ylabel('Objective ','FontSize',16);
return

function y=eval_fun(x)
y=((6*x-2).^2).*sin(12*x-4);
return
