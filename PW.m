function W= PW(N)
W=zeros(1,N);
n=log2(N);
%PW���й��켫���룬��������޹�
for i=0:N-1
    w(i+1,:)=bitget(i,1:n);
end
for iii=1:N
    for j=1:n
    W(iii)=w(iii,j)*(2.^((1/4)*(j-1)))+W(iii);
    end
end