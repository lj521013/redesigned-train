%****************ģ���˻��㷨��С��********************
function [Xval1,Yval1]=SimulateAnnealing1
clear all;
close all;
clc;
tic;
tem=30; 
step = 0.02; 
Xm= 1;%x���ֵ
Ym = 1;%y���ֵ
err = 1e-8; 
mlen = 10000;% ����ɷ�������
sj = 0.95; %˥������
accp = 0.0; 
rnd =rand;
disp('*******ģ���˻��㷨�����н������*******:');
disp(sprintf('��ʼ�¶�: %d',tem));
disp(sprintf('˥������: %d',sj));
disp(sprintf('����ɷ�������: %d',mlen));
disp(sprintf('��������: %d',step));
disp(sprintf('��������Ϊ: �������Ž�֮��С��%d',err));
% ���ѡ�� ��ֵ�趨
Xp = -Xm * rand; %���ó�ʼֵ
Yp = -Ym * rand;
Xval = Xp; 
Yval = Yp;
Xp = -Xm * rand; 
Yp = -Ym * rand;
Xval1 = Xp; 
Yval1 = Yp;
mm=abs( f( Xval1,Yval1)-f (Xval,Yval));%�������Ž�֮��
%ģ���˻��㷨
while mm > err
    tem=sj*tem;
    accp = 0.0;  
    for i=0:mlen:1 %����mlen��    
        p=0;
        while p==0
            Xnext = Xp + step*Xm*(rand-0.5);%���ѡ��һ��
            NextY = Yp+ step*Ym*(rand-0.5);          
            if p==(~(Xnext >= -Xm && Xnext <= Xm && NextY>= -Ym && NextY <= Ym))
                p=1;
            end
        end
        if (f(Xval1,Yval1) >f(Xnext,NextY))%�Ƿ�ȫ�����Ž�
            Xval =Xval1;
            Yval = Yval1;          
            Xval1=Xnext;%�µ����Ž�
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
disp(sprintf('����ֵ��Ϊ: [%.4f,%.4f]',Xval1,Yval1));
disp(sprintf('����ֵΪ: %d',val));
runtime=toc;
disp(sprintf('����ʱ��Ϊ: %d',runtime));
end
%****************�Ż�����********************
function  value=f(x,y)
          value=-cos(2.*x+1).*cos(2.*y+1)-4.*cos(3.*x+2).*cos(3.*y+2)-9.*cos(4.*x+3).*cos(4.*y+3)-16.*cos(5.*x+4).*cos(5.*y+4)-25.*cos(6.*x+5).*cos(6.*y+5);
end