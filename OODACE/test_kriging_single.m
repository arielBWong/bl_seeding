function test_kriging
clear all;close all;clc;
addpath E:\BING\OODACE
startup;

x_trg = linspace(0, 1, 10); %[0,0.5,1]';
y_trg = eval_fun(x_trg);


x_trge = [0.05, 0.55, 0.95];
y_trge = eval_fun(x_trge);

k = oodacefit(x_trg,y_trg);

x_tst = [linspace(0,1,100)]';
y_tst=eval_fun(x_tst);
[y_prd,~] = k.predict(x_tst);

plot(x_tst,y_tst,'r-');hold on;
plot(x_tst,y_prd,'b-');
plot(x_trg,y_trg,'ro','Markersize',15); hold on;


function y= eval_fun(x)
y=((6*x-2).^2).*sin(12*x-4);
return
