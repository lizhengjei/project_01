function out=erfc(x)
syms t;
y=exp(-t^2);
out=int(y,t,x,inf);%int����������֣�x���ޣ�inf���ޣ���������y
out=2*out/sqrt(pi);
clear t
end