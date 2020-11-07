clear all
close all 
%initialisation
t=1:1000;
deltat=1;

xpoint = 5; %m/s
ypoint = -7 ; %m/s
v=[xpoint, ypoint];% vecteur vitesse du projectile

sigma=10*pi/180; % ecart type
epsilon=sigma*randn(1,length(t));


Zsb=zeros(1,length(t));
z=zeros(1,length(t));
Ksi=zeros(2,length(t));
Ksik=zeros(2,length(t));

%P=zeros(2,2,length(t));
Hk=zeros(2,length(t));
%K=zeros(2,length(t));

p=zeros(2,2,length(t));
X=zeros(1,length(t));
Y=zeros(1,length(t));

F=[1 0 ; 0 1];

X(1)=1000;
Y(1)=1000;

Zsb(1)=atan2(X(1),Y(1)); % signal sans bruit 
z(1)=Zsb(1)+epsilon(1);  % signal bruité 

rg=1400; % distance intuitive (la distance que l'on pense entre nous et l'objet)

Ksik(:,1)=rg*[sin(z(1)) cos(z(1))];

p=10^6*[1 0 ; 0 1];

for i=2:length(t)
    X(i)=X(i-1)+deltat*v(1);
    Y(i)=Y(i-1)+deltat*v(2);
    P=F*p*F';
    
    Zsb(i)=atan2(X(i),Y(i));
    z(i)=Zsb(i)+epsilon(i);

    Ksi(:,i)=F*Ksik(:,i-1)+deltat*v';
    
    
    Hk=1/(Ksi(1,i)^2+Ksi(2,i)^2)*[Ksi(2,i) -Ksi(1,i)];
    K=P*Hk'*inv(Hk*P*Hk'+sigma^2); 

    Ksik(:,i)=Ksi(:,i)+K*(z(i)-atan2(Ksi(1,i),Ksi(2,i)));
    p=P-K*Hk*P;
end

%  Affichage des résulats
figure(1)
plot(Ksik(1,:),Ksik(2,:), 'b')
hold on
plot(X,Y, 'g')
legend('prédiction','trajectoire réel')
xlabel('x (en m)')
ylabel('y (en m)')

figure(2)

plot(t,Zsb, 'g','linew',5)
hold on 
plot(t,z, 'r')
legend('angle sans bruit','angle bruité')
xlabel('t (en s)')
ylabel('angle (en rad)')
