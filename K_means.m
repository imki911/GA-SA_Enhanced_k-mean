%2018年6月6日
%标准k均值聚类
clear all
clc
 [x1 x2 x3 x4 type]=textread('bezdekIris.txt','%f,%f,%f,%f,%s');
% [x1 x2 x3 x4 type]=textread('moni.txt','%f,%f,%f,%f,%s');
X=[x1 x2 x3 x4];
sampleTotal=size(x1,1);
X(1:50,5)=1;
X(51:100,5)=2;
X(101:150,5)=3;

%乱序排列输入样本
randNum=rand(sampleTotal,1);
[mm,n]=sort(randNum);% n是整数序号
Xdisorder=X(n,:);

scaleMin=min(X(:,1:4));
scaleMax=max(X(:,1:4));

% center1=scaleMin+(scaleMax-scaleMin).*rand(1,4);
% center2=scaleMin+(scaleMax-scaleMin).*rand(1,4);
% center3=scaleMin+(scaleMax-scaleMin).*rand(1,4);

center1=X(ceil(150*rand(1)),1:4);
center2=X(ceil(150*rand(1)),1:4);
center3=X(ceil(150*rand(1)),1:4);

center11=1.1*center1;
center22=1.1*center2;
center33=1.1*center3;
count=0

while(sum(center11~=center1) || sum(center22~=center2) || sum(center33~=center3) )
count=count+1;
class1=[];
class2=[];
class3=[];

for i=1:150
    dist1=norm((center1-X(i,1:4)));
    dist2=norm((center2-X(i,1:4)));
    dist3=norm((center3-X(i,1:4)));
    if dist1<dist2 && dist1< dist3
        class1=[class1;X(i,1:5)];
    elseif dist2<dist1 && dist2< dist3
        class2=[class2;X(i,1:5)];
    else 
        class3=[class3;X(i,1:5)]; 
    end
    

    
end
    center11=center1;
    center1=mean(class1(:,1:4),1);
    
    center22=center2;
    center2=mean(class2(:,1:4),1);
    
    center33=center3;
    center3=mean(class3(:,1:4),1);

end %end while
%%
%计算fitness
a=0.1;
b=0.5;
Sigma_in=0;
J_b=0;
mean_all=mean(X(:,1:4));
J_b=sum((center1-mean_all).^2)+ sum((center2-mean_all).^2 ) ...
    + sum((center3-mean_all).^2 );

Sigma_bw=sum((center1-center2).^2)+sum((center1-center3).^2)...
    +sum((center2-center3).^2);

Sigma_in=Sigma_in+sum(sum((class1(:,1:4)-repmat(center1,size(class1,1),1)).^2,1 ),2 );
Sigma_in=Sigma_in+sum(sum((class2(:,1:4)-repmat(center2,size(class2,1),1)).^2 ,1),2);
Sigma_in=Sigma_in+sum(sum((class3(:,1:4)-repmat(center3,size(class3,1),1)).^2 ,1 ),2);
Sigma_in=Sigma_in/150;
% fitness=Sigma_bw/(a+b*Sigma_in)
fitness=J_b/Sigma_in
count;
