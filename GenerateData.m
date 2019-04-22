clear all
clc 
close all
u1=[1.2 3.4 5.2 6.4];
u2=[5.0 3.1 5.2 4.2];
u3=[2.6 5.2 1.3 4.5 ];
sigma1=[2 0 0 0;
        0 3 0 0; 
        0 0 1 0;
        0 0 0 4 ];
sigma2=[3 0 0 0;
        0 1 0 1;
        0 0 10 3; 
        0 1 3 8];
sigma3=[4 0 0 0; 
        0 2 0 0; 
        0 0 12 9; 
        0 0 9 11];
X1=mvnrnd(u1,sigma1,50);
X2=mvnrnd(u2,sigma2,50);
X3=mvnrnd(u3,sigma3,50);
X=[X1;X2;X3];

%存储仿真结果数据
fid=fopen('moni.txt','w');%建立文件
%循环写入数据
for i=1:150
    fprintf(fid,'%.3f,%.3f,%.3f,%.3f,xxx\n',X(i,1),X(i,2),X(i,3),X(i,4));%保存小数点后8位，Ipluse_data为结果数据  
end
fclose(fid);