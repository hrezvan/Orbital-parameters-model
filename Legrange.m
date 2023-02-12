function [a,e,i,omega,wp,f]=Legrange(A,t)

T=A(:,1);

kk=zeros(7,1);
s=zeros(8,1);
P_n=zeros(6,1);

for i=2:7
    P_t=A(:,i);
    for j=1:8
        n=[1:8];
        n(n==j)=[];
        for k=1:size(n,2)
            kk(k) = (t-T(n(k)))/(T(j)-T(n(k)));
        end
        s(j)=P_t(j)*prod(kk);
    end
    P_n(i-1) =sum(s);
end

r= sqrt(P_n(1)^2+P_n(2)^2+P_n(3)^2);
Vs=sqrt(P_n(4)^2+P_n(5)^2+P_n(6)^2);
miu= 398600 *(1000)^3 ;
a=(miu*r)/(2*miu-r*(Vs^2));

sin_eta=(P_n(1)*P_n(4)+P_n(2)*P_n(5)+P_n(3)*P_n(6))/(r*Vs);

e = ((Vs^2-miu/r)*[P_n(1) P_n(2) P_n(3)]-dot([P_n(1) P_n(2) P_n(3)],[P_n(4) P_n(5) P_n(6)])*[P_n(4) P_n(5) P_n(6)])/miu;
f=acos(dot(e,[P_n(1) P_n(2) P_n(3)])/(norm(e)*r));
h=cross(P_n(1:3),P_n(4:6));
i= acos(h(3)/norm(h));

n = cross([0 0 1],h);
wp = acos(dot(n,e)/(norm(n)*norm(e)));
if e(3)<0
   wp = 2*pi-wp;
end

omega = acos(n(1)/norm(n));
if n(2)<0
   omega = 2*pi-omega;
end
e = norm(e);

end
