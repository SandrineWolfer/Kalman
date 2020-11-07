close all
clear all
%Initialisation 
t=1:500; 
Omega=1; 

F=[2*cos(Omega) -1 ; 1 0];
Hk=[1 0];

epsilon=2;
psi=0;

z=zeros(1,length(t));    
x=zeros(2,length(t));

x(:,1)=[pi/3;0];

a=1;
Zsb(t)=a*cos(Omega*t+psi); % signal sans bruit 
z(t)=Zsb(t)+epsilon*randn(1,length(t));   % signal bruité 

y(1)=Hk*x(:,1)+epsilon*randn;

p=10^6*[1 0 ; 0 1];%p0 sachant 0

% Boucle filtre de Kalman 
for k=2:length(t)
    Xkk=F*x(:,k-1);
    P=F*p*F'
    K=P*Hk'*inv(Hk*P*Hk'+epsilon^2);
    x(:,k)=Xkk+K*(z(k)-Hk*Xkk);
    p=P-K*Hk*P;
end

%  Affichage des résulats
figure
plot(z,'r')
hold on
plot(x(1,:),'b')
hold on
plot(Zsb,'g')
legend('bruité','prédiction','non bruité')
xlabel('temps')
ylabel('amplitude')


