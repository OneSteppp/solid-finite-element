m=20;%将（0,1）分为m个区间
h=1/m;%步长
u=zeros(m+1,1);%待求解量
B=zeros(m+1,1);
A=zeros(m+1,m+1);%系数矩阵
x=0:1/m:1;%x轴
syms Xi;%计算标准化，坐标变换
for i=1:m
    f11=(1+(x(i)+h*Xi)^2)/h+(sin(x(i)+h*Xi)^2)*((1-Xi)^2)*h;
    a11=int(f11,Xi,0,1);%积分
    f12=-(1+(x(i)+h*Xi)^2)/h+(sin(x(i)+h*Xi)^2)*Xi*(1-Xi)*h;
    a12=int(f12,Xi,0,1);
    f22=(1+(x(i)+h*Xi)^2)/h+(sin(x(i)+h*Xi)^2)*(Xi^2)*h;
    a22=int(f22,Xi,0,1);
    %
    A(i,i)=A(i,i)+a11;
    A(i+1,i+1)=A(i+1,i+1)+a22;
    A(i+1,i)=A(i+1,i)+a12;
    A(i,i+1)=A(i,i+1)+a12;
    %
    f1=exp(x(i)+h*Xi)*(1-Xi)*h;
    f2=exp(x(i)+h*Xi)*Xi*h;
    b1=int(f1,Xi,0,1);
    b2=int(f2,Xi,0,1);
    %
    B(i)=B(i)+b1;
    B(i+1)=B(i+1)+b2;
end
%处理边界条件
A(1,1)=10e20;
A(m+1,m+1)=10e20;%u(0)=0,u(1)=0
%追赶法
%LU分解
a=zeros(m+1,1);
b=zeros(m+1,1);
c=zeros(m+1,1);
a(1)=A(1,1);
b(1)=A(1,2)/A(1,1);
for i=2:m+1
    c(i)=A(i,i-1);%L矩阵的下三角部分
end
for i=2:m
    b(i)=A(i,i+1)/(A(i,i)-A(i,i-1)*b(i-1));%U矩阵的上三角部分
end
for i=2:m+1
    a(i)=A(i,i)-A(i,i-1)*b(i-1);%L矩阵的对角线部分
end
%求解Ly=B
y=zeros(m+1,1);
y(1)=B(1)/a(1);
for i=2:m+1
    y(i)=(B(i)-c(i)*y(i-1))/a(i);
end
%求解Uu=y
u(m+1)=y(m+1);
for i=1:m
    u(m+1-i)=y(m+1-i)-b(m+1-i)*u(m+2-i);
end
plot(x,u);
xlabel('x');
ylabel('u');


