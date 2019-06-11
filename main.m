%数据准备
EE=2.1*10^5; %弹性模量，N/mm2
uu=0.3; %泊松比
hh=10; %厚度，mm
cx=40; %节点x坐标基准
cy=20; %节点y坐标基准
cj=zeros(55,2); %创建节点坐标矩阵
m=0; %计数器
for i=55:-5:5 %给cj矩阵赋值
    k=0;
    for j=i:-1:i-4
        cj(j,1)=m*cx;
        cj(j,2)=k*cy;
        k=k+1;
    end
    m=m+1;
end
m=0;k=0;i=0;j=0; %计数器归零
cd=zeros(80,3); %创建单元节点编号矩阵
m=1;
for i=1:8:80 %给cd矩阵赋值
    k=m;
    for j=i:i+3
        cd(j,1)=k;
        cd(j+4,1)=k;
        cd(j,2)=k+6;
        cd(j+4,2)=k+5;
        cd(j,3)=k+1;
        cd(j+4,3)=k+6;
        k=k+1;
    end
    m=m+5;
end
m=0;k=0;i=0;j=0; %计数器归零
pp(110)=0; %载荷矩阵，N
pp(2)=-200; %载荷矩阵赋值
pp(102)=-200;
for i=12:10:92 %载荷矩阵中间各点赋值
    pp(i)=-400;
end
i=0; %计数器归零
cysw=zeros(110,1); %约束位移值
cysj=[51 52 53 53 55]; %约束位移节点
 
%刚度矩阵计算
kk=zeros(110,110);
l=0;m=0;n=0; %计数器
for i=1:80 %遍历各单元计算单元刚度矩阵，直接赋值进入整体刚度矩阵
    k0=zeros(110,110); %创建临时刚度矩阵
    clear a b c;
    l=cd(i,1);%赋值单元编号
    m=cd(i,2);
    n=cd(i,3);
    delta=0.5*(cj(l,1)*cj(m,2)+cj(m,1)*cj(n,2)+cj(n,1)*cj(l,2))-0.5*(cj(m,1)*cj(l,2)+cj(n,1)*cj(m,2)+cj(l,1)*cj(n,2)); %计算参数
    a(l)=(cj(m,1)*cj(n,2)-cj(m,2)*cj(n,1));
    b(l)=-(cj(n,2)-cj(m,2));
    c(l)=(cj(n,1)-cj(m,1));
    a(m)=-(cj(l,1)*cj(n,2)-cj(n,1)*cj(l,2));
    b(m)=(cj(n,2)-cj(l,2));
    c(m)=-(cj(n,1)-cj(l,1));
    a(n)=(cj(l,1)*cj(m,2)-cj(m,1)*cj(l,2));
    b(n)=-(cj(m,2)-cj(l,2));
    c(n)=(cj(m,1)-cj(l,1));
    for j=1:3 %r,s=l,m,n
        for k=j:3
            if j==1
            r=l;
            else if j==2
              r=m;
                 else r=n;
                end
            end
            if k==1
            s=l;
            else if k==2
                 s=m;
                 else s=n;
                end
            end
            k0(2*r-1,2*s-1)=(EE*hh)*(b(r)*b(s)+((1-uu)*c(r)*c(s))/2)/(4*(1-uu^2)*delta);
            k0(2*s-1,2*r-1)=k0(2*r-1,2*s-1);
            k0(2*r,2*s-1)=(EE*hh)*(uu*c(r)*b(s)+((1-uu)*b(r)*c(s))/2)/(4*(1-uu^2)*delta);
            k0(2*s-1,2*r)=k0(2*r,2*s-1);
            k0(2*r-1,2*s)=(EE*hh)*(uu*b(r)*c(s)+((1-uu)*c(r)*b(s))/2)/(4*(1-uu^2)*delta);
            k0(2*s,2*r-1)=k0(2*r-1,2*s);
            k0(2*r,2*s)=(EE*hh)*(c(r)*c(s)+((1-uu)*b(r)*b(s))/2)/(4*(1-uu^2)*delta);
            k0(2*s,2*r)=k0(2*r,2*s);
        end
    end
    kk=kk+k0; %将临时刚阵赋值给整体刚阵
end
 
%根部位移约束处理
for i=1:length(cysj)
    pp(cysj(i)*2)=0;
    pp(cysj(i)*2-1)=0;
end
for i=101:110
    kk(i,i)=10e15;
end
 
%节点位移求解
[RA,RB,n,dd]=liezhu(kk,pp');
 
%绘图
size_huitu=20; %绘图变形放大倍数控制：20倍
for i=1:55 %对节点位移值进行格式调整
dd_huitu(i,1)=dd(2*i-1);
dd_huitu(i,2)=dd(2*i);
end
cj_huitu=cj+size_huitu*dd_huitu; %将节点变形前坐标和变形位移叠加得新节点坐标
x=cj_huitu(:,1); %节点变形后x坐标
y=cj_huitu(:,2); %节点变形后y坐标
figure(1) %绘变形前后对比图
plot(x,y,'.b') %绘制变形后图形
for i=1:80 %绘制三角形
m1=x(cd(i,1));
m2=y(cd(i,1));
n1=x(cd(i,2));
n2=y(cd(i,2));
p1=x(cd(i,3));
p2=y(cd(i,3));
line([m1,n1],[m2,n2],'color','b')
line([n1,p1],[n2,p2],'color','b')
line([p1,m1],[p2,m2],'color','b')
end
hold on
plot(cj(:,1),cj(:,2),'.r') %绘制变形前图形
for i=1:80 %绘制三角形
m1=cj(cd(i,1),1);
m2=cj(cd(i,1),2);
n1=cj(cd(i,2),1);
n2=cj(cd(i,2),2);
p1=cj(cd(i,3),1);
p2=cj(cd(i,3),2);
line([m1,n1],[m2,n2],'color','r')
line([n1,p1],[n2,p2],'color','r')
line([p1,m1],[p2,m2],'color','r')
end
xlim([-1,500])
ylim([-100,100])
title('悬臂梁变形前后对比图（变形放大20倍）')
hold off
clear x y  s r p2 p1 n2 n1 m2 m1 m n l
clear k0 i j k delta c RA RB%清除过程变量
