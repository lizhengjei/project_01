% function index = getIndex(start, last, M)
% index = zeros(1,M);
% sum = start + last;
% index(M) = last;
% for i = 1:M-1
%     if i==1
%         index(i) = floor((sum-1)/(M-1)) + start;
%     else
%        index(i) = floor((sum-1)/(M-1)) + index(i-1);
%     end
% end
% end
% Lp=4;
% PM=[0.4,1,4,7,3,0.8,10,20];
% active_index=[2,4,5,7,8];
%  [~,index]=sort(PM(active_index));
% PM_min_index = active_index(index(1:Lp))
% A=5;
% ismember(A,PM_min_index)
Lp=4;
len = 2;
activepath = [1 0 1 1 1 0 1 0];
active_index = [3 4];
PM=[0.9 1 2 4 2 1 9 1 ];
min_index=find(activepath==1)%�ڴ���·���н�PM����
                    [ia, ~] = setdiff(min_index, active_index)%�ڴ��·���н�ͨ����żУ���·��ȥ��
%                     res = min_index(sort(ia))
                     [~,min_index_]=sort(PM(ia))%��ʣ��·������PM��С��Lp-len��·������
               min=ia(min_index_(1:Lp-len))
               real_index = [active_index min]
 