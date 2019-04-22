%2018年6月7日
%遗传模拟退火算法
%考虑使用自适应的变异参数，适应度函数设置有问题

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
% 输入样本绘图
% for i=1:150
%     plot(Xdisorder(i,1:4));
%     hold on
% end

%----参数设置--------
m=50;%种群规模
pc=0.8;%交叉概率
pm=0.2;%变异概率
popu=zeros(m,sampleTotal);%种群的染色体表示
T=420;
%----参数设置--------

%初始化种群
for i=1:m
    popu(i,:)=ceil(3*rand(1,sampleTotal));
end

popuNew=popu;

count=0;
fitnessBestRec=[];%记录最好的类内方差
fitnessMinRec = [];
fitnessEveRec=[];
fitnessMin=0;
fitnessEverage=100;
fitnessMax=inf;
%while(count<200 || (fitnessMax-fitnessEverage)>0.01)
while(count<150 || fitnessBestRec(count-105)~= fitnessBestRec(count-1))
count=count+1;
if T >10 
    T=T-1;%温度降低
end
%选择
if  size(popu,1)>m %若上一次进化之后种群规模大于预设
    fitnessRecord=zeros(size(popu,1),1);%存储总体类间方差
    %fitnessRecord=zeros(size(popu,1),1);%存储适应度
    for i=1:size(popu,1)
        fitnessRecord(i)=fitness(popu(i,:),Xdisorder(:,1:4));
    end
    [bw_m, bw_n]=sort(fitnessRecord,'descend');
    popuNew10pencent=popu((bw_n(1:ceil(0.1*m) )),:);%保留前10%的最优个体基因
    
    %剩下的个体轮盘赌
    popuremainRoulette=popu( (bw_n(1+ceil(0.1*m):end )),:);%剩余的用于轮盘赌的个体
    remainIndividual=size(popuremainRoulette,1);%剩余的个体总数
    remainToBeGeneratedIndividual=m-ceil(0.1*m); %轮盘赌中要选出的个体数量
    addedFlag=zeros(remainIndividual,1);%标记该个体已经被选中过了
    fitnessSum=sum(fitnessRecord(bw_n(ceil(0.1*m)+1:end)) );
    fitnessRouletteProb= fitnessRecord(bw_n(ceil(0.1*m)+1:end) )/fitnessSum; %剩余个体对应的轮盘赌概率
    fitnessRouletteProbMax=max(fitnessRouletteProb);%最大适应度个体的轮盘赌概率
    
    fitnessMax=max( fitnessRecord);%(ceil(0.1*m)+1:end) );
    fitnessMin=min( fitnessRecord);%(ceil(0.1*m)+1:end) );
    fitnessEverage=mean(fitnessRecord);%(ceil(0.1*m)+1:end));
    
    
    
    
    popuNewRoulette=zeros(remainToBeGeneratedIndividual,size(popu,2));%初始化存储轮盘赌产生个体的矩阵
    RouletteGenerated=0;%已经产生的个体计数
    loopIndex=0;
    while (RouletteGenerated<remainToBeGeneratedIndividual)
         index=mod(loopIndex,remainIndividual)+1; %对剩余个体的下标循环
         if 3*fitnessRouletteProbMax*rand(1)<fitnessRouletteProb(index) && addedFlag(index)==0 
             RouletteGenerated=RouletteGenerated+1;
             popuNewRoulette(RouletteGenerated,:)=popuremainRoulette(index,:);
             addedFlag(index)=1;
         end
         loopIndex=loopIndex+1;
    end
    
    popuBest=popu((bw_n(1)),:);
    fitnessRecord(bw_n(1));%显示最优个体的适应度
    fitnessBestRec=[fitnessBestRec;fitnessRecord(bw_n(1))];
    fitnessMinRec=[fitnessMinRec;fitnessMin];
    fitnessEveRec=[fitnessEveRec;fitnessEverage];
    popuNew=[popuNew10pencent;popuNewRoulette];
    popu=popuNew;
end
%交叉


for i=1:round(m/2)
    if rand(1)<pc
        fatherA=popu(  ceil(m*rand(1)),: );
        fatherB=popu(  ceil(m*rand(1)),: );        
        pos=ceil(60*rand(1))+2;
        child1=[fatherA(1:pos) fatherB(pos+1:150)];
        child2=[fatherB(1:pos) fatherA(pos+1:150)];
        
%         fatherA=child1;
%         fatherB=child2;
%         pos=ceil(70*rand(1))+50;
%         child1=[fatherA(1:pos) fatherB(pos+1:150)];
%         child2=[fatherB(1:pos) fatherA(pos+1:150)];
        popu=[popu;child1;child2];
    end
end
%变异
for i=2:size(popu,1)
    if rand(1)<pm
        fitnessVariOld=fitness(popu(i,:),Xdisorder(:,1:4));
        temp1=popu(i,:);
        temp2=popu(i,:);
        pos=ceil(sampleTotal*rand(1));
        temp1(pos)=mod(temp1(pos),3)+1;
        temp2(pos)=mod(temp2(pos)+1,3)+1;

        fitnessVariNew1 =fitness(temp1,Xdisorder(:,1:4));
        fitnessVariNew2 =fitness(temp2,Xdisorder(:,1:4));
        
        
        if fitnessVariNew1 > fitnessVariOld
            popu=[popu;temp1];
        elseif rand(1)< exp((fitnessVariNew1-fitnessMin)/(fitnessEverage-fitnessMin)/T)
                popu(i,:)=temp1;
        end
        
         if fitnessVariNew1 > fitnessVariOld
            popu=[popu;temp2];
        elseif rand(1)< exp((fitnessVariNew1-fitnessMin)/(fitnessEverage-fitnessMin)/T)
                popu=[popu;temp2];
         end       
        
        
%         popu(i,:)=temp;
        
    end
end

end  %end while

%%
%效果检验
sampleLabelOrderRecover=zeros(150,1);
for i=1:150
   % sampleLabelOrderRecover(n(i),:)=Xdisorder(i,:);
   sampleLabelOrderRecover(n(i))=popuBest(i);
end
figure
bar(sampleLabelOrderRecover)
xlabel('迭代次数')
ylabel('类别')
figure
plot(fitnessBestRec,'r') ; hold on
plot(fitnessMinRec,'g') ;hold on
plot(fitnessEveRec,'b');
legend('Max','Min','Everage')
xlabel('迭代次数')
ylabel('适应度')
legend('Max','Min','Everage')