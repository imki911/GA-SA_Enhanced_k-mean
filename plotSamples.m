clear all
clc
 [x1 x2 x3 x4 type]=textread('bezdekIris.txt','%f,%f,%f,%f,%s');
%[x1 x2 x3 x4 type]=textread('moni.txt','%f,%f,%f,%f,%s');
X=[x1 x2 x3 x4];
sampleTotal=size(x1,1);
X(1:50,5)=1;
X(51:100,5)=2;
X(101:150,5)=3;
plot(X(1:50,1:4)','r');
hold on
plot(X(51:100,1:4)','g');
hold on
plot(X(101:150,1:4)','b');
legend('c1','c2','c3')
title('数据集2'),xlabel('坐标轴')