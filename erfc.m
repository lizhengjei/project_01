function out=erfc(x)
syms t;
y=exp(-t^2);
out=int(y,t,x,inf);%int函数计算积分，x下限，inf上限，被积函数y
out=2*out/sqrt(pi);
clear t
end