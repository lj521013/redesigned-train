%****************�Ŵ��㷨��С��********************
function[X]=hanshuga2 
clear all;
close all;
clc; 
tic;
pc=0.9;%�������
pm=0.05;%�������
n=200;%Ⱦɫ�����
iter=200; %�Ŵ�����
u=rst(n,44,0,0,0); %����100*44��Ⱦɫ�����
u1=u(:,1:22);
u2=u(:,23:44);
[N,L]=size(u1);
Xm=-1;xmax=1;%��Сֵ�����ֵ
Ym=-1;ymax=1;
x0=-1:0.02:1;
y0=-1:0.02:1;
[X,Y]=meshgrid(x0,y0);  
f1=-cos(2.*X+1).*cos(2.*Y+1)-4.*cos(3.*X+2).*cos(3.*Y+2)-9.*cos(4.*X+3).*cos(4.*Y+3)-16.*cos(5.*X+4).*cos(5.*Y+4)-25.*cos(6.*X+5).*cos(6.*Y+5);
val=1;
aver=[ ];%����������ֵ
it=1;%�Ŵ�������
xval=[ ];yval=[ ];%�����������ֵ
C=[ ];
x=jm(u1,Xm,xmax);
y=jm(u2,Ym,ymax);
f='-cos(2.*x+1).*cos(2.*y+1)-4.*cos(3.*x+2).*cos(3.*y+2)-9.*cos(4.*x+3).*cos(4.*y+3)-16.*cos(5.*x+4).*cos(5.*y+4)-25.*cos(6.*x+5).*cos(6.*y+5)';
fv=eval(f);%���㺯��ֵ
%��ͼ
set(0,'defaultfigurecolor','w')
figure(1);%Ⱦɫ���ʼλ��ͼ1
mesh(X,Y,f1,'LineWidth',1);grid on;hold on;plot3(x,y,fv,'b *');
xlabel('x');ylabel('y');zlabel('f');hold off;
figure(2);%Ⱦɫ���ʼλ��ͼ2
contour(X,Y,f1,'LineWidth',1);grid on;hold on;plot3(x,y,fv,'b *');
xlabel('x');ylabel('y');zlabel('f');hold off;
pb=0.5;
pw=0.1;
pr=1-(pb+pw); %���Ƹ���   
nb=round(pb * N);
nw=round(pw * N);
nr=round(pr * N);
if(nb+nw+nr)~=N, %��֤���ƺ�Ⱦɫ����������
    dif=N-(nb+nw+nr);
    nr=nr+dif;
end;
%���ɳ�ʼȺ��
while it <=iter     
    [rw,ind]=sort(fv);%����
    ind=fliplr(ind);%�������ҷ�ת
    vt1=[u1(ind(1:nb),:);u1(ind(end-nw+1:end),:);u1(2:nr+1,:)];%���Ʋ���
    vt2=[u2(ind(1:nb),:);u2(ind(end-nw+1:end),:);u2(2:nr+1,:)];
    %�������
    C(:,1)=rand(N,1)<=pc;%����50*1����
    C(:,2)=round(22. * rand(N,1));
    I=find(C(:,1)==1);
    IP=[I,C(I,2)];
    %vtemp�ǳ�ʼȺ�弯�ϣ��������,���ɽ����
    for i=1:size(IP,1),
        u1(IP(i,1),:)=[vt1(IP(i,1),1:IP(i,2)) vt1(1,IP(i,2)+1:end)];
        u2(IP(i,1),:)=[vt2(IP(i,1),1:IP(i,2)) vt2(1,IP(i,2)+1:end)];
    end
    %�������
    m1=rand(N,L)<=pm; %����50*22����
    m2=rand(N,L)<=pm;
    m1(1,:)=zeros(1,L); %��һ�е�Ⱦɫ�岻��
    m2(1,:)=zeros(1,L); 
        u1=u1-2.*(u1.*m1)+m1;
        u2=u2-2.*(u2.*m2)+m2;
        x=jm(u1,Xm,xmax);%����
        y=jm(u2,Ym,ymax);       
        fv=-eval(f);
        [val,indb]=max(fv);%�洢��it��������ֵ
        u1(1,:)=u1(indb,:);%�洢����Ⱦɫ��
        u2(1,:)=u2(indb,:);
        media=mean(fv);%ƽ��ֵ
        xval=[xval val]; 
        yval=[yval val];
        aver=[aver media];  
        it=it+1;           
end;
disp(sprintf('*******�Ŵ��㷨�����н������*******'));
disp(sprintf('�����: %d',iter));
disp(sprintf('�˿�: %d',N));
disp(sprintf('�������: %.3f',pc));
disp(sprintf('ͻ�����: %.3f',pm));
disp(sprintf('����ֵ��Ϊ: [%.4f,%.4f]',x(indb),y(indb)));
disp(sprintf('����ֵΪ: %d',-val));
%��ͼ
figure(3);%Ⱦɫ������λ��ͼ1
mesh(X,Y,f1,'LineWidth',1);grid on;hold on;plot3(x,y,fv,'b*');
xlabel(('x'));ylabel('y');zlabel('f'); 
figure(4);%Ⱦɫ������λ��ͼ2
contour(X,Y,f1,'LineWidth',1);grid on;hold on;plot3(x,y,fv,'b*');
xlabel(('x'));ylabel('y');zlabel('f'); 
figure(5);%����ֵ�仯ͼ��
plot(-xval,'r *','LineWidth',1);xlabel('Generations');ylabel('value');
figure(6)%ƽ��ֵ�仯ͼ��
plot(-aver,'b','LineWidth',1);xlabel('Generations');ylabel('average value');hold off;
runtime=toc
disp(sprintf('����ʱ��Ϊ: %d',runtime));

function x=jm(v1,xmin,xmax); %decode����
v1=fliplr(v1);%��������ҷ�ת
s=size(v1);%
aux=0:1:21;aux=ones(s(1),1)*aux;
x1=sum((v1.*2.^aux)');
x=xmin+(xmax-xmin)*x1./4194303;



%**********����Ⱦɫ�����(����*����)************
function[m,n]=rst(n1,s1,n2,s2,bip)
if nargin==2  
    n2=n1;s2=s1;bip=1; 
  elseif nargin==4, 
    bip=1;  
end;
m=2.*rand(n1,s1)-1; %rand(n1,s1)��ʾ����n1 X s1��ֵΪ0-1��������� 
if bip==1,
    m=hardlims(m);
else;
    m=hardlim(m);
end;
n=2.*rand(n2,s2)-1;
if bip==1,
    n=hardlims(n);
else;
    n=hardlim(n);
end;



