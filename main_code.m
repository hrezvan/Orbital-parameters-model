clc;
clear;
close all;

%% input data
xls1=xlsread('coordinate_group2.xlsx');
xls2=xlsread('ephemeris_group2.xlsx');
ground=xls1(:,1:4) ;
image(:,1)=xls1(:,1);
image(:,2:3)=xls1(:,5:6);
ephemeris=xls2(:,2:8);
Nr=6000;
f=1.084; 
PS=13*10^-6;
beta=-20.84*pi/180;     %% Cross-track angle
t_cl=377.373;

%% Parameters
[a,e,inc,w0,wp,f0]=Legrange(ephemeris,t_cl);
w1=1.991*10^-7;

X=zeros(15,1);
X(1)=f0;
X(3)=w0;
X(4)=w1;
X(5)=a;
X(6)=inc;

par=[beta,wp,e,f];
%% points

x=(image(:,3)-(Nr/2)).*PS; 
y=(image(:,2)-(Nr/2)).*PS;

GCP_Points=[26 12 24 23 38 36 35 27 6 7 2 33 20 18 10 15 17];

points = [x y image(:,3) ground(:,2) ground(:,3) ground(:,4)];
id=1:size(points,1);
id(GCP_Points)=[];
ICP_Points=id;
check_point=points(ICP_Points,:);
control_point=points(GCP_Points,:);                                

%% structure
% str = [1 1 1 1 1 1 1 0 0 1 0 0 1 0 0];
% str = [1 1 1 1 1 1 0 0 0 0 0 0 0 0 0];
% str = [1 1 1 1 1 1 1 1 0 1 1 0 1 1 0];
str = [1 1 1 1 1 1 0 1 0 0 1 0 0 1 0];
% str = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];

%% Case two
B_quasi=-eye(sum(str==1));
opm2=[X X];
dX=1;
count=2;
while count<100
 F2=opm2(str==1,count-1)-opm2(str==1,count);
 for i=1:size(control_point,1)
    [Be(2*i-1:2*i,:),F1(2*i-1:2*i,1)]=B_F(control_point(i,:),par,X,str);
 end
 B=[Be;B_quasi];
 F=[F1;-F2];

 P=eye(length(F));
 
 dX=inv(B'*P*B)*B'*P*F;
 X(str==1)=X(str==1)+dX;
 opm2(:,count+1)=X;
 count=count+1;
end

for i=1:size(check_point,1)
[XY_c(i,1:2) Dxy(i,1:2)]=RMSE(check_point(i,:),par,X);
end
RMSE_xy=sqrt(sum((Dxy(:,1).^2+Dxy(:,2).^2))/length(check_point(:,1)));
disp ('RMSE image (meter) = ');disp((RMSE_xy));
RMSE_XY=sqrt(sum(((check_point(:,4)-XY_c(:,1)).^2+(check_point(:,5)-XY_c(:,2)).^2))/length(check_point(:,1)));
disp ('RMSE ground (meter) = ');disp((RMSE_XY));

figure
subplot(1,2,1)
plot(points(:,4), points(:,5),'^r');
hold on
plot(XY_c(:,1),XY_c(:,2),'ob');
quiver(XY_c(:,1), XY_c(:,2), (check_point(:,4)-XY_c(:,1)), (check_point(:,5)-XY_c(:,2)));
t=title('Ground points');

subplot(1,2,2)
plot(control_point(:,1),control_point(:,2),'^r');
hold on
plot(check_point(:,1),check_point(:,2),'ob');
quiver(check_point(:,1),check_point(:,2), Dxy(:,1), Dxy(:,2));
t=title('Image points');
