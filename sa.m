%****************模拟退火算法求极小点********************
function [Xval1,Yval1]=SimulateAnnealing1
clear all;
close all;
clc;
tic;
tem=30; 
step = 0.02; 
Xm= 1;%x最大值
Ym = 1;%y最大值
err = 1e-8; 
mlen = 10000;% 马尔可夫链长度
sj = 0.95; %衰减参数
accp = 0.0; 
rnd =rand;
disp('*******模拟退火算法的运行结果如下*******:');
disp(sprintf('初始温度: %d',tem));
disp(sprintf('衰减参数: %d',sj));
disp(sprintf('马尔可夫链长度: %d',mlen));
disp(sprintf('步长因子: %d',step));
disp(sprintf('结束条件为: 两次最优解之差小于%d',err));
% 随机选点 初值设定
Xp = -Xm * rand; %设置初始值
Yp = -Ym * rand;
Xval = Xp; 
Yval = Yp;
Xp = -Xm * rand; 
Yp = -Ym * rand;
Xval1 = Xp; 
Yval1 = Yp;
mm=abs( f( Xval1,Yval1)-f (Xval,Yval));%两次最优解之差
%模拟退火算法
while mm > err
    tem=sj*tem;
    accp = 0.0;  
    for i=0:mlen:1 %迭代mlen次    
        p=0;
        while p==0
            Xnext = Xp + step*Xm*(rand-0.5);%随机选下一点
            NextY = Yp+ step*Ym*(rand-0.5);          
            if p==(~(Xnext >= -Xm && Xnext <= Xm && NextY>= -Ym && NextY <= Ym))
                p=1;
            end
        end
        if (f(Xval1,Yval1) >f(Xnext,NextY))%是否全局最优解
            Xval =Xval1;
            Yval = Yval1;          
            Xval1=Xnext;%新的最优解
            Yval1=NextY;           
        end             
        if( f(Xp,Yp)- f(Xnext,NextY) >0 ) 
            Xp=Xnext;
            Yp=NextY;
            accp=accp+1;
        else
            changer = -1 * ( f(Xnext,NextY) -f(Xp,Yp) ) / tem ;
            rnd=rand;
            p1=exp(changer);
            double (p1);
            if p1 >rand 
                Xp=Xnext;
                Yp=NextY;
                accp=accp+1;
            end
        end
    end
    mm=abs( f( Xval1,Yval1)-f (Xval,Yval));
    val=f(Xval1, Yval1);
end
disp(sprintf('最优值点为: [%.4f,%.4f]',Xval1,Yval1));
disp(sprintf('最优值为: %d',val));
runtime=toc;
disp(sprintf('运行时间为: %d',runtime));
end
%****************优化函数********************
function  value=f(x,y)
          value=-cos(2.*x+1).*cos(2.*y+1)-4.*cos(3.*x+2).*cos(3.*y+2)-9.*cos(4.*x+3).*cos(4.*y+3)-16.*cos(5.*x+4).*cos(5.*y+4)-25.*cos(6.*x+5).*cos(6.*y+5);
end