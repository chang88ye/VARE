function obj=subjective_te(weight,ind,idealpoint,D)
%�������һ������������һ����ֵ��
%weight��20*2��
%ֻ�ʺ���һ������
s=size(weight,1);  %niche
indsize=size(ind,1); %������niche����Ҳ������һ��

weight((weight==0))=0.00001;

if indsize==s 
   % size(ind(:,D+1:Dim))
  %  idealpoint
     %size(idealpoint(ones(1,indsize),:))
%ind(:,D+1,Dim)-idealpoint(ones(1,indsize),:)
    part2=abs(ind(:,D+1:end)-idealpoint(ones(1,indsize),:));
    obj=max(weight.*part2,[],2);%%Ҫ������������һ��������
elseif indsize==1   %��ֻ��һ������������
    part2=abs(ind(D+1:end)-idealpoint);
    
    obj=max(weight.*part2(ones(1,s),:),[],2);
else
    error('�����ǲ���');
end