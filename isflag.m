function flag = isflag
x=rand(1,1);
if x>=0.5
    flag =1;
else 
    flag = 0;
% x=rand(1,100);
% yes=0;
% no=0;
% for i = 1: 100
%     if x(i)<0.5
%         yes=yes+1;
%     else
%         no=no+1;
%     end
% end
% if yes >= no
%     flag = 1;
% else
%     flag = 0;
% end

end