function fitness = fitness( individual,sample )
%����һ��Ⱦɫ�壬ԭʼ������������Ӧ�ȣ��������ڷ��
%   Detailed explanation goes here
global mean_all;
global class1Center;
global class2Center
global class3Center
global J_b;
global class1;
global class2;
global class3;

class=3;%������
Sigma_in=0;%��ʼ�����ڷ����С����
Sigma_bw=0; %��ʼ����䷽��
J_b=0;%���ƽ����׼��(���)
a=0.1;
b=0.5;
sampleTotal=150;%����������Ⱦɫ�峤�ȣ�
class1=[];
class2=[];
class3=[];

    for j=1:sampleTotal
        if individual(j)==1
            class1=[class1;sample(j,:)];
        elseif individual(j)==2
            class2=[class2;sample(j,:)];
        elseif individual(j)==3
            class3=[class3;sample(j,:)];
        end %end if
    end %end for
    
    
    mean_all=mean(sample,1);
    
    %�������ڷ���
    if ~isempty(class1)  
    class1Center=mean(class1(:,1:4),1);;
    Sigma_in=Sigma_in+sum(sum((class1-repmat(class1Center,size(class1,1),1)).^2,1 ) ,2);
    J_b=J_b+ sum((class1Center-mean_all).^2 );
    end
    if ~isempty(class2) 
    class2Center=mean(class2(:,1:4),1);
    Sigma_in=Sigma_in+sum(sum((class2-repmat(class2Center,size(class2,1),1)).^2 ));
    J_b=J_b+ sum((class2Center-mean_all).^2 );
    end
    if ~isempty(class3) 
    class3Center=mean(class3(:,1:4),1);
    Sigma_in=Sigma_in+sum(sum((class3-repmat(class3Center,size(class3,1),1)).^2 ));
    J_b=J_b+ sum((class3Center-mean_all).^2  );
    end
    
%    J_b=sum(sum((class1Center-mean_all).^2,1))+sum(sum((class2Center-mean_all).^2,1))...
%        +sum(sum((class3Center-mean_all).^2,1));
    
    %����䷽��
    if ~isempty(class1) && ~isempty(class2)
        Sigma_bw=Sigma_bw+sum((class1Center-class2Center).^2);
    end
    if ~isempty(class1) && ~isempty(class3)
        Sigma_bw=Sigma_bw+sum((class1Center-class3Center).^2);
    end
    if ~isempty(class2) && ~isempty(class3)
        Sigma_bw=Sigma_bw+sum((class2Center-class3Center).^2);
    end
    
        
  
    
%     Sigma_bw=sum((class1Center-class2Center).^2)+sum((class1Center-class3Center).^2)...
%         +sum((class2Center-class3Center).^2);
    
%     Sigma_in=Sigma_in+sum(sum((class1-repmat(class1Center,size(class1,1),1)).^2 ) );
%     Sigma_in=Sigma_in+sum(sum((class2-repmat(class2Center,size(class2,1),1)).^2 ));
%     Sigma_in=Sigma_in+sum(sum((class3-repmat(class3Center,size(class3,1),1)).^2 ));

 Sigma_in=Sigma_in/150;
%     fitness=Sigma_bw/(a+b*Sigma_in);
   fitness= J_b/Sigma_in;
   % fitness=1000/Sigma_in;
    %Sigma_bw;
    
    
end

