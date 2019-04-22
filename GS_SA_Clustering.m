%2018��6��7��
%�Ŵ�ģ���˻��㷨
%����ʹ������Ӧ�ı����������Ӧ�Ⱥ�������������

clear all
clc
  [x1 x2 x3 x4 type]=textread('bezdekIris.txt','%f,%f,%f,%f,%s');
% [x1 x2 x3 x4 type]=textread('moni.txt','%f,%f,%f,%f,%s');
X=[x1 x2 x3 x4];
sampleTotal=size(x1,1);
X(1:50,5)=1;
X(51:100,5)=2;
X(101:150,5)=3;

%����������������
randNum=rand(sampleTotal,1);
[mm,n]=sort(randNum);% n���������
Xdisorder=X(n,:);
% ����������ͼ
% for i=1:150
%     plot(Xdisorder(i,1:4));
%     hold on
% end

%----��������--------
m=50;%��Ⱥ��ģ
pc=0.8;%�������
pm=0.2;%�������
popu=zeros(m,sampleTotal);%��Ⱥ��Ⱦɫ���ʾ
T=420;
%----��������--------

%��ʼ����Ⱥ
for i=1:m
    popu(i,:)=ceil(3*rand(1,sampleTotal));
end

popuNew=popu;

count=0;
fitnessBestRec=[];%��¼��õ����ڷ���
fitnessMinRec = [];
fitnessEveRec=[];
fitnessMin=0;
fitnessEverage=100;
fitnessMax=inf;
%while(count<200 || (fitnessMax-fitnessEverage)>0.01)
while(count<150 || fitnessBestRec(count-105)~= fitnessBestRec(count-1))
count=count+1;
if T >10 
    T=T-1;%�¶Ƚ���
end
%ѡ��
if  size(popu,1)>m %����һ�ν���֮����Ⱥ��ģ����Ԥ��
    fitnessRecord=zeros(size(popu,1),1);%�洢������䷽��
    %fitnessRecord=zeros(size(popu,1),1);%�洢��Ӧ��
    for i=1:size(popu,1)
        fitnessRecord(i)=fitness(popu(i,:),Xdisorder(:,1:4));
    end
    [bw_m, bw_n]=sort(fitnessRecord,'descend');
    popuNew10pencent=popu((bw_n(1:ceil(0.1*m) )),:);%����ǰ10%�����Ÿ������
    
    %ʣ�µĸ������̶�
    popuremainRoulette=popu( (bw_n(1+ceil(0.1*m):end )),:);%ʣ����������̶ĵĸ���
    remainIndividual=size(popuremainRoulette,1);%ʣ��ĸ�������
    remainToBeGeneratedIndividual=m-ceil(0.1*m); %���̶���Ҫѡ���ĸ�������
    addedFlag=zeros(remainIndividual,1);%��Ǹø����Ѿ���ѡ�й���
    fitnessSum=sum(fitnessRecord(bw_n(ceil(0.1*m)+1:end)) );
    fitnessRouletteProb= fitnessRecord(bw_n(ceil(0.1*m)+1:end) )/fitnessSum; %ʣ������Ӧ�����̶ĸ���
    fitnessRouletteProbMax=max(fitnessRouletteProb);%�����Ӧ�ȸ�������̶ĸ���
    
    fitnessMax=max( fitnessRecord);%(ceil(0.1*m)+1:end) );
    fitnessMin=min( fitnessRecord);%(ceil(0.1*m)+1:end) );
    fitnessEverage=mean(fitnessRecord);%(ceil(0.1*m)+1:end));
    
    
    
    
    popuNewRoulette=zeros(remainToBeGeneratedIndividual,size(popu,2));%��ʼ���洢���̶Ĳ�������ľ���
    RouletteGenerated=0;%�Ѿ������ĸ������
    loopIndex=0;
    while (RouletteGenerated<remainToBeGeneratedIndividual)
         index=mod(loopIndex,remainIndividual)+1; %��ʣ�������±�ѭ��
         if 3*fitnessRouletteProbMax*rand(1)<fitnessRouletteProb(index) && addedFlag(index)==0 
             RouletteGenerated=RouletteGenerated+1;
             popuNewRoulette(RouletteGenerated,:)=popuremainRoulette(index,:);
             addedFlag(index)=1;
         end
         loopIndex=loopIndex+1;
    end
    
    popuBest=popu((bw_n(1)),:);
    fitnessRecord(bw_n(1));%��ʾ���Ÿ������Ӧ��
    fitnessBestRec=[fitnessBestRec;fitnessRecord(bw_n(1))];
    fitnessMinRec=[fitnessMinRec;fitnessMin];
    fitnessEveRec=[fitnessEveRec;fitnessEverage];
    popuNew=[popuNew10pencent;popuNewRoulette];
    popu=popuNew;
end
%����


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
%����
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
%Ч������
sampleLabelOrderRecover=zeros(150,1);
for i=1:150
   % sampleLabelOrderRecover(n(i),:)=Xdisorder(i,:);
   sampleLabelOrderRecover(n(i))=popuBest(i);
end
figure
bar(sampleLabelOrderRecover)
xlabel('��������')
ylabel('���')
figure
plot(fitnessBestRec,'r') ; hold on
plot(fitnessMinRec,'g') ;hold on
plot(fitnessEveRec,'b');
legend('Max','Min','Everage')
xlabel('��������')
ylabel('��Ӧ��')
legend('Max','Min','Everage')