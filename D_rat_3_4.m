%% Izhikevich   4:1 AB  电耦合规则连接   10层每层100个   放电率随耦合度变化
clear;clc;
tic
m = 10;        % 10层
n = 5;         % 每层5个
n1=n*0.8;          % A number
n2=n*0.2;          % B number
w=cell(m,n);   % 初始化胞组
a = [ones(1,n1)*0.04  ones(1,n2)*0.15];
b = [ones(1,n1)*0.2   ones(1,n2)*0.17];
c = [ones(1,n1)*(-65)  ones(1,n2)*(-65)];
d = [ones(1,n1)*10  ones(1,n2)*2];
dt = 0.5;
tspan = 0:dt:100;    
amp = 0.5;e=1;
for jj=0:0.1:1
    g = ones(1,n)* jj;
    for t=tspan                          
       if(t>0)                        
          I = ones(1,n)*100;
       else
          I=zeros(1,n); 
       end 
       V=cell(m,n); 
       u=cell(m,n); 
       VF=cell(m,n);  
         for i=1:m                      % m 行  3
             for k = 1:n                  % n 列 5
                 for j=1:length(tspan)-1
                     V{i,k}(1,1)=-65;
                     u{i,k}(1,1)=0.2*(-65);
                     if i==1   
                         V{i,k}(1,j+1)=V{i,k}(1,j)+dt*(0.04*V{i,k}(1,j)^2+5*V{i,k}(1,j)+140-u{i,k}(1,j))+I(k); % first lay stimulaton
                         u{i,k}(1,j+1)=u{i,k}(1,j)+dt*a(k)*(b(k)*V{i ,k}(1,j)-u{i,k}(1,j));
                     else
                         V{i,k}(1,j+1)=V{i,k}(1,j)+dt*(0.04*V{i,k}(1,j)^2+5*V{i,k}(1,j)+140-u{i,k}(1,j))-g(k)*(V{i,k}(1,j)-V{i-1,k}(1,j)); % other lay couper first lay
                         u{i,k}(1,j+1)=u{i,k}(1,j)+dt*a(k)*(b(k)*V{i ,k}(1,j)-u{i,k}(1,j));
                     end
                      if V{i,k}(1,j+1)>30                                                                %膜电位大于30mV时，辅助复位机制
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
    for i=1:m                                 % m 行   
             for k = 1:n                  % n 列   
                  sp_n =V{i,k}(1,:); 
                  sp_m = findpeaks(sp_n,'minpeakdistance',1,'minpeakheight',-40);
                  sp_m1 = length(sp_m); 
                  spk_rat(i,k) = sp_m1/max(tspan);
             end
%               plot(spk_rat(i,1),'-'); hold on
    end
%     for jj=1:1:m
%         A_rat1_i(e)=spk_rat(jj,n);
%     end
  
A_rat1(e)=spk_rat(1,n1);  %兴奋性神经元
A_rat2(e)=spk_rat(2,n1);
A_rat3(e)=spk_rat(3,n1);
A_rat4(e)=spk_rat(4,n1);
A_rat5(e)=spk_rat(5,n1);
A_rat6(e)=spk_rat(6,n1);
A_rat7(e)=spk_rat(7,n1);
A_rat8(e)=spk_rat(8,n1);
A_rat9(e)=spk_rat(9,n1);
A_rat10(e)=spk_rat(10,n1);
%   A_rat4(e)=spk_rat(1,4);A_rat5(e)=spk_rat(2,4);A_rat6(e)=spk_rat(3,4);
A_rat11(e)=spk_rat(1,n);  %抑制性神经元
A_rat12(e)=spk_rat(2,n);
A_rat13(e)=spk_rat(3,n);
A_rat14(e)=spk_rat(4,n);
A_rat15(e)=spk_rat(5,n);
A_rat16(e)=spk_rat(6,n);
A_rat17(e)=spk_rat(7,n);
A_rat18(e)=spk_rat(8,n);
A_rat19(e)=spk_rat(9,n);
A_rat20(e)=spk_rat(10,n);
  e=e+1;
end     
D=0:0.1:1;
figure(2)
% plot(D,A_rat4,'-diamonk','Markersize',8,'LineWidth',1.5);hold on    %  兴奋性神经元
% plot(D,A_rat5,'-.or','Markersize',8,'LineWidth',1.5);hold on
% plot(D,A_rat6,'-->b','Markersize',8,'LineWidth',1.5);
plot(D,A_rat11,'-diamonk','Markersize',8,'LineWidth',1.5);hold on   %  抑制性神经元
plot(D,A_rat12,'-.or','Markersize',8,'LineWidth',1.5);hold on
plot(D,A_rat13,'-->b','Markersize',8,'LineWidth',1.5);hold on
plot(D,A_rat14,'-.og','Markersize',8,'LineWidth',1.5);hold on
plot(D,A_rat15,'-.oc','Markersize',8,'LineWidth',1.5);hold on
plot(D,A_rat16,'-.om','Markersize',8,'LineWidth',1.5);hold on
plot(D,A_rat17,'-.oy','Markersize',8,'LineWidth',1.5);hold on
plot(D,A_rat18,'-.ok','Markersize',8,'LineWidth',1.5);hold on
plot(D,A_rat19,'-.om','Markersize',8,'LineWidth',1.5);hold on
plot(D,A_rat20,'-.ob','Markersize',8,'LineWidth',1.5);
legend({'lay1','lay2','lay3','lay4','lay5','lay6','lay7','lay8','lay9','lay10'});  % 加 图 例   
% legend({'lay1','lay2','lay3'});
xlabel('电耦合度','FontSize',30,'fontname','times new roman');
ylabel('频率/kHz','FontSize',30,'fontname','times new roman');
set(gca,'FontSize',30);set(gca,'ylim',[0 1]); 
grid on
set(gca,'GridLineStyle','--','GridColor','k', 'GridAlpha',0.5);
set(gca,'Fontname', 'Times New Roman','FontSize',30);
set(gcf,'color',[1,1,1]);  

figure(1)
plot(D,A_rat1,'-diamonk','Markersize',8,'LineWidth',1.5);hold on   %  抑制性神经元
plot(D,A_rat2,'-.or','Markersize',8,'LineWidth',1.5);hold on
plot(D,A_rat3,'-->b','Markersize',8,'LineWidth',1.5);hold on
plot(D,A_rat4,'-.og','Markersize',8,'LineWidth',1.5);hold on
plot(D,A_rat5,'-.oc','Markersize',8,'LineWidth',1.5);hold on
plot(D,A_rat6,'-.om','Markersize',8,'LineWidth',1.5);hold on
plot(D,A_rat7,'-.oy','Markersize',8,'LineWidth',1.5);hold on
plot(D,A_rat8,'-.ok','Markersize',8,'LineWidth',1.5);hold on
plot(D,A_rat9,'-.om','Markersize',8,'LineWidth',1.5);hold on
plot(D,A_rat10,'-.ob','Markersize',8,'LineWidth',1.5);
legend({'lay1','lay2','lay3','lay4','lay5','lay6','lay7','lay8','lay9','lay10'});  % 加 图 例                    
xlabel('电耦合度','FontSize',30,'fontname','times new roman');                           % X 轴标注设置字体大小
ylabel('频率/kHz','FontSize',30,'fontname','times new roman');                       % Y 轴标注设置字体大小
set(gca,'ylim',[0 1]);                               % Y 轴设置显示范围
set(gca,'FontSize',30);                              % 坐标轴刻度设置大小
grid on
set(gca,'GridLineStyle','--','GridColor','k', 'GridAlpha',0.5);
set(gca,'Fontname', 'Times New Roman','FontSize',30);
set(gcf,'color',[1,1,1]);  
%% 画A/B神经元每层放电波形
% for j=n+2:n+3
%     figure(j)
%      y=ones(1,length(tspan));
%     for i=0:1:m-1
%         plot3(tspan,y+i,VF{i+1,j-3}(1,:));hold on
%     end  
%     xlabel('Time/ms');ylabel('lay');zlabel('Voltage/mV');set(gca,'YTick',1:1:m);
%     grid on
%     set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
%     j=j+1;
% end
%%
toc
disp(['runtime: ',num2str(toc)]);