function c = calcu(p)
% c=cell(1,M);
temp=[];
% for i = 1:M
if length(p)<=4
    temp=p;
else
q=1;
    for j =1:length(p)
        flag=isflag;
        if flag==1
            temp(q)=p(j);
            q=q+1;
        end
    end
end
c=temp;
end