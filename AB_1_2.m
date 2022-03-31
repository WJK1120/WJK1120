%% Izhikevich   8:2 AB  ����Ϲ�������   10��ÿ��100��  ��ʾ������ԪĤ��λ����
clear;clc;
tic
m = 10;        % 10��
n = 5;         % ÿ��5��
n1=n*0.8;          % A number
n2=n*0.2;          % B number
w=cell(m,n);   % ��ʼ������
a = [ones(1,n1)*0.04  ones(1,n2)*0.15];
b = [ones(1,n1)*0.2   ones(1,n2)*0.17];
c = [ones(1,n1)*(-65)  ones(1,n2)*(-65)];
d = [ones(1,n1)*10  ones(1,n2)*2];

g = ones(1,n)* 0.2;                 % �����ϵ�� 0.2
% Vsyn = 2;sigma = 10; theta = -0.35;

% dt = 0.25;
dt = 0.5;
tspan = 0:dt:100;    
amp = 0.5;

for t=tspan                          
       if(t>0)                        
          I = ones(1,n)*100;                          % ���̼���һ����Ԫ
%         x=0:0.1:10; I = gaussmf(x,[2 5]);I = I*100; % ģ�ⳬ���ֲ��̼�          
       else
          I=zeros(1,n); 
       end 
       V=cell(m,n); 
       u=cell(m,n); 
       VF=cell(m,n);  
         for i=1:m                      % m ��  3
             for k = 1:n                  % n �� 5
                 for j=1:length(tspan)-1
                     V{i,k}(1,1)=-65;
                     u{i,k}(1,1)=0.2*(-65);
                     if i==1   
                         V{i,k}(1,j+1)=V{i,k}(1,j)+dt*(0.04*V{i,k}(1,j)^2+5*V{i,k}(1,j)+140-u{i,k}(1,j))+I(k);  % first lay stimulaton
                         u{i,k}(1,j+1)=u{i,k}(1,j)+dt*a(k)*(b(k)*V{i ,k}(1,j)-u{i,k}(1,j));
                     else
                         V{i,k}(1,j+1)=V{i,k}(1,j)+dt*(0.04*V{i,k}(1,j)^2+5*V{i,k}(1,j)+140-u{i,k}(1,j))-g(k)*(V{i,k}(1,j)-V{i-1,k}(1,j)); % other lay couper first lay
                         u{i,k}(1,j+1)=u{i,k}(1,j)+dt*a(k)*(b(k)*V{i ,k}(1,j)-u{i,k}(1,j));
                     end
                      if V{i,k}(1,j+1)>30                                                                %Ĥ��λ����30mVʱ��������λ����
                              VF{i,k}(1,j+1)=30;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                              V{i,k}(1,j+1)=c(k);
                              u{i,k}(1,j+1)=u{i,k}(1,j+1)+d(k);                 
                      else              
                              VF{i,k}(1,j+1)=V{i,k}(1,j);
                      end
                 end 
             end
         end
end

sp_n = zeros(1,length(tspan)-1);
spk_rat = zeros(m,n);  
%  /////////////////  rat  //////////////
% for i=1:m                                 % m ��   
%              for k = 1:n                  % n ��   
%                   sp_n =V{i,k}(1,:); 
%                   sp_m = findpeaks(sp_n,'minpeakdistance',1,'minpeakheight',-40);
%                   sp_m1 = length(sp_m); 
%                   spk_rat(i,k) = sp_m1/max(tspan);
%              end
% %               plot(spk_rat(i,1),'-'); hold on
% end
%% �ֿ���ʾÿ����ԪĤ��λ
% for j=1:1:n
%   figure (j) 
%     for i=1:1:m
%         subplot(m,1,i); plot(tspan,VF{i,j}(1,:)) ;
%     end
% end
%% ������ԪĤ��λ
figure(n+1)
plot(spk_rat(:,1),'-'); xlabel('lay');ylabel('rat');ylim([0 1.4]);set(gca,'XTick',1:1:m);
%% ��ά��ʾһ���˷�����Ԫ��һ����������Ԫ
for j=n+2:n+3
    figure(j)
     y=ones(1,length(tspan));
    for i=0:1:m-1
        plot3(tspan,y+i,VF{i+1,j-3}(1,:));hold on         % ��ά��ʾÿ��Ĥ��λ
    end
    xlabel('ʱ��/ms','FontSize',20,'fontname','times new roman');
    ylabel('����','FontSize',20,'fontname','times new roman');
    zlabel('Ĥ��λ/mV','FontSize',20,'fontname','times new roman');
    set(gca,'YTick',1:1:m);
    grid on
    set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
    set(gca,'FontSize',20);
    set(gca,'zlim',[-100 50]); 
    j=j+1;
    view(40,70);set(gcf,'color',[1,1,1]);
end
toc
disp(['runtime: ',num2str(toc)]);