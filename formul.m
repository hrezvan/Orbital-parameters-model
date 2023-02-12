clc
clear all

syms W0 W1 f0 f1 a in e beta wp w0 w1 w2 Phi0 Phi1 Phi2 k0 k1 k2 X Y Z  t c x y x0 y0 r f

fi=f0+f1*t;
W=W0+W1*t;

r=(a*(1-(e^2)))/(1+e*cos(fi));

ws=w0+w1*t+w2*(t.^2);
ps=Phi0+Phi1*t+Phi2*(t.^2);
ks=k0+k1*t+k2*(t.^2);

R_coor1=[1           0         0
         0 cos(in - pi/2) sin(in - pi/2)
         0 -sin(in - pi/2) cos(in - pi/2)];
  
R_coor2=[cos(pi/2 - (fi+wp)) 0 -sin(pi/2 - (fi+wp))
              0          1        0
         sin(pi/2 - (fi+wp)) 0 cos(pi/2 - (fi+wp))];

R_coor3=[cos(pi - W) sin(pi - W) 0
        -sin(pi - W) cos(pi - W) 0
             0          0        1];
      
R1=[1           0               0
    0 cos(-(in - pi/2)) sin(-(in - pi/2))
    0 -sin(-(in - pi/2)) cos(-(in - pi/2))];

R2=[cos(-(pi/2 - (fi+wp))) 0 -sin(-(pi/2 - (fi+wp)))
            0              1         0
    sin(-(pi/2 - (fi+wp))) 0 cos(-(pi/2 - (fi+wp)))];
  
R3=[cos(-(pi - W)) sin(-(pi - W)) 0
   -sin(-(pi - W)) cos(-(pi - W)) 0
           0             0        1];

Rw=[1     0     0
    0 cos(ws) sin(ws)
    0 -sin(ws) cos(ws)];

Rp=[cos(ps) 0 -sin(ps)
      0     1     0
    sin(ps) 0 cos(ps)];

Rk=[cos(ks) sin(ks) 0
   -sin(ks) cos(ks) 0
       0      0     1];

R1_beta=[1    0        0
         0 cos(beta) sin(beta)
         0 -sin(beta) cos(beta)];

R_coor=R_coor3*R_coor1*R_coor2;
coor=R_coor*[0;0;r];
coor=simplify(coor);
Xs=coor(1,1);Ys=coor(2,1);Zs=coor(3,1);

R=R2*R1*R3;
R=simplify(R);

R_wpk=Rk*Rp*Rw;

DX=[X-Xs
    Y-Ys
    Z-Zs];

C=R1_beta*R_wpk*R*DX;

Fx=f*(C(1)/C(3));
Fy=y+f*(C(2)/C(3));

B=jacobian([Fx Fy],[f0 f1 W0 W1 a in w0 w1 w2 Phi0 Phi1 Phi2 k0 k1 k2]);

dx=[x
    y
   -f];

R_2=R_coor3'*R_coor1'*R_coor2'*Rw'*Rp'*Rk'*R1_beta';
C2=R_2*dx;

X_C=(C2(1)/C2(3))*(Z-Zs)+Xs;
Y_C=(C2(2)/C2(3))*(Z-Zs)+Ys;


