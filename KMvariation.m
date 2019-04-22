function [ chromVariated ] = KMvariation( chrom,Sample)
%k均值变异算子
%根据当前染色体求距离中心；将距离聚类中心最近的基因位改写为该类
%   
class1=[];
class2=[];
class3=[];
for i=1:150
    class=chrom(i);
    if class==1
    class1=[class1; Sample(i,1:4)];
    elseif class==2
    class2=[class2;Sample(i,1:4)];
    else
    class3=[class3;Sample(i,1:4)];
    end
end
%     center1=mean(class1);
%     center2=mean(class2);
%     center3=mean(class3);
%     size(Sample(:,1:4) ) 
%    size( repmat(center1,150,1) )


    distance2class1=inf(150,1);
    distance2class2=inf(150,1);
    distance2class3=inf(150,1);
    
    if ~isempty(class1)  
        center1=mean(class1,1);
        
        distance2class1=sum((Sample(:,1:4)-repmat(center1,150,1)).^2,2);
    end
    if ~isempty(class2)  
        center2=mean(class2,1);
        distance2class2=sum((Sample(:,1:4)-repmat(center2,150,1)).^2,2);
    end
    if ~isempty(class3)  
        center3=mean(class3,1);
        distance2class3=sum((Sample(:,1:4)-repmat(center3,150,1)).^2,2);
    end
    
    chromVariated=chrom;
    
    %最小距离分类
    for i=1:150
        if distance2class1(i)<distance2class2(i) && distance2class1(i)<distance2class3(i)
            chromVariated(i)=1;
        elseif distance2class2(i)<distance2class1(i) && distance2class2(i)<distance2class3(i)
            chromVariated(i)=2;
        else
            chromVariated(i)=3;
        end
        
    end
    
%     [m,n]=sort(distance2class1);
%     chromVariated(n(1))=1;
%     
%     [m,n]=sort(distance2class2);
%     chromVariated(n(1))=2; 
%     
%     [m,n]=sort(distance2class3);
%     chromVariated(n(1))=3;
        
    
end

