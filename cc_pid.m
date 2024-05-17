clear all;
close all;
syms x;
%% PID 
%%
Y=zeros(10,1);
I_N = eye(5);
I_m = eye(1);
Ee_PID = zeros(5,1);
e1_PID = zeros(5,1);
%  [ yy ]=woa_cc_pid(2);
%  kp=yy(1,1)*eye(5)
%  ki=yy(1,2)*eye(5)
%  kd=yy(1,3)*eye(5)
%% 
for k=1:1:121
time(k)=(k-1)/20;
%% 定义系统参数
t=0.35;
A1 = t* [ 1     1
            -1.55  0.65 ];
A2 = t* [ 1     1
            -2  1.1 ];
B = t*[ 0;
          1 ];  
A1_d = t* [   0.07  0.02;
             0.09  0.051 ];  
A2_d = t* [    0.05  0.078;
             0.05  0.017 ];  
    
C = [ 0.01  1 ];

L = [ 2  -1  -1  0  0;
      -1  2  -1  0  0;
     -1 0  3   -1  -1;
      -1 -1  -1   4   -1;
     -1  0 -1  -1   3];
M = [ 1 0 0 0 0;
      0 1 0 0 0;
      0 0 1 0 0;
      0 0 0 1 0;
      0 0 0 0 1 ];
H = L + M;
    
    % 参考信号
    if k<15
       y_R1(:,k) =zeros(5,1);
    elseif k>=15 & k<=30
        y_R1(:,k) = ( (k-15)/10 )*100*ones(5,1);
    else
       y_R1(:,k) = 1.5*100*ones(5,1);
    end
    
     if k<15
       y_R(:,k) =zeros(5,1);
    elseif k>=15 & k<=30
        y_R(:,k) = ( (k-15)/10 )*100*ones(5,1);
    else
       y_R(:,k) = 1.5*100*ones(5,1);
     end
    
     if k+1<15
       y_R(:,k+1) =zeros(5,1);
    elseif k+1 >=15 & k+1 <=30
        y_R(:,k+1) = ( (k+1-15)/10 )*100*ones(5,1);
    else
       y_R(:,k+1) = 1.5*100*ones(5,1);
     end
    
     if k+2<15
       y_R(:,k+2) =zeros(5,1);
    elseif k+2 >=15 & k+2 <=30
        y_R(:,k+2) = ( (k+2-15)/10 )*100*ones(5,1);
    else
       y_R(:,k+2) = 1.5*100*ones(5,1);
     end
    
     if k+3<15
       y_R(:,k+3) =zeros(5,1);
    elseif k+3 >=15 & k+3 <=30
        y_R(:,k+3) = ( (k+3-15)/10 )*100*ones(5,1);
    else
       y_R(:,k+3) = 1.5*100*ones(5,1);
     end
    
    Y1(:,k)=[Y(1,k) Y(3,k) Y(5,k) Y(7,k) Y(9,k)];
    Y2(:,k)=[Y(2,k) Y(4,k) Y(6,k) Y(8,k) Y(10,k)];
%     Y(:,k)=[ Y(1,k);Y(2,k);Y(3,k);Y(4,k);Y(5,k);Y(6,k);Y(7,k);Y(8,k);Y(9,k);Y(10,k) ];
    Tr1(:,k)=kron(eye(5),C)*Y(:,k);
    %% 定义模糊代号
     h1(:,k) =  ( 20*(   sin( (Y(2,k)+Y(4,k)+Y(6,k)+Y(8,k)+Y(10,k))/5  )   )^2  ) / (exp(2)+exp(5));
     h2(:,k) = 1-h1(:,k);
     A_pi_0 = h1(k)*A1 + h2(k)*A2;
     Y(:,k+1) = kron(eye(5),A_pi_0)*Y(:,k);
     
     h1(k+1) = ( 20*(   sin( (Y(2,k+1)+Y(4,k+1)+Y(6,k+1)+Y(8,k+1)+Y(10,k+1))/5  )   )^2  ) / (exp(2)+exp(5));
     h2(k+1) = 1-h1(k+1);
     A_pi_1 = h1(k+1)*A1 + h2(k+1)*A2;
     Y(:,k+2) = kron(eye(5),A_pi_1)*Y(:,k+1);
     
     h1(k+2) = ( 20*(   sin( (Y(2,k+2)+Y(4,k+2)+Y(6,k+2)+Y(8,k+2)+Y(10,k+2))/5  )   )^2  ) / (exp(2)+exp(5));
     h2(k+2) = 1-h1(k+2);
     A_pi_2 = h1(k+2)*A1 + h2(k+2)*A2;
     
     Ad_pi_0 = h1(k)*A1_d + h2(k)*A2_d;
     Ad_pi_1 = h1(k+1)*A1_d + h2(k+1)*A2_d;
     Ad_pi_2 = h1(k+2)*A1_d + h2(k+2)*A2_d;
     
     %e_PID(:,k) = [ e_1(k); e_2(k); e_3(k) ];
g(:,k) = 0.004*[abs(Y(1,k)+1) - abs(Y(1,k)-1); abs(Y(3,k)+1) - abs(Y(3,k)-1); abs(Y(5,k)+1) - abs(Y(5,k)-1); abs(Y(7,k)+1) - abs(Y(7,k)-1); abs(Y(9,k)+1) - abs(Y(9,k)-1) ];
    
g(:,k+1) = 0.004*[abs(Y(1,k+1)+1)-abs(Y(1,k+1)-1); abs(Y(3,k+1)+1)-abs(Y(3,k+1)-1); abs(Y(5,k+1)+1)-abs(Y(5,k+1)-1); abs(Y(7,k+1)+1)-abs(Y(7,k+1)-1); abs(Y(9,k+1)+1)-abs(Y(9,k+1)-1) ];
    
g(:,k+2) = 0.004*[abs(Y(1,k+2)+1)-abs(Y(1,k+2)-1); abs(Y(3,k+2)+1)-abs(Y(3,k+2)-1); abs(Y(5,k+2)+1)-abs(Y(5,k+2)-1); abs(Y(7,k+2)+1)-abs(Y(7,k+2)-1); abs(Y(9,k+2)+1)-abs(Y(9,k+2)-1) ];

    
    W(:,k) = t*[ 0;  -0.01*cos(Y(1,k));  0;  -0.01*cos(Y(3,k));  0;  -0.01*cos(Y(5,k));  0;  -0.01*cos(Y(7,k)); 0;  -0.01*cos(Y(9,k)) ];
    % 误差 e
    e_PID(:,k)= -kron(H,I_m) * (kron(I_N,C)*Y(:,k) - y_R(:,k));
 
      kp = [ 0.15         0         0         0         0;
             0         0.1232         0         0         0;
             0              0    0.145         0         0;
             0              0         0    0.1232         0;
             0              0         0         0    0.1232 ];

     ki = [ 1.25         0         0         0         0;
                 0    1.1616         0         0         0;
                 0         0    1.18         0         0;
                 0         0         0    1.1616         0;
                 0         0         0         0    1.1616 ];

     kd = [ 0.75         0         0         0         0;
                 0    0.6668         0         0         0;
                 0         0    0.666         0         0;
                 0         0         0    0.6668         0;
                 0         0         0         0    0.6668 ];

     %% I[u_k]
     u(:,k) = kp*e_PID(:,k) + ki*Ee_PID + kd*(e_PID(:,k) - e1_PID );
     Ee_PID = Ee_PID + e_PID(:,k);
     e1_PID = e_PID(:,k);
     
     if k >3
     Y(:,k+1) =  h1(:,k) * (  kron(I_N,A1)*Y(:,k) + kron(I_N,A1_d)*Y(:,k-3) + kron(I_N,B)* ( u(:,k)...
                                            + 0.004*[abs(Y(1,k)+1) - abs(Y(1,k)-1); abs(Y(3,k)+1) - abs(Y(3,k)-1); abs(Y(5,k)+1) - abs(Y(5,k)-1); ...
                                                     abs(Y(7,k)+1) - abs(Y(7,k)-1); abs(Y(9,k)+1) - abs(Y(9,k)-1) ] )  + ...
                            t*[ 0;  -0.01*cos(Y(1,k));  0;  -0.01*cos(Y(3,k));  0;  -0.01*cos(Y(5,k));  0;  -0.01*cos(Y(7,k)); 0;  -0.01*cos(Y(9,k)) ] ) +...
                 h2(:,k) * (  kron(I_N,A2)*Y(:,k) + kron(I_N,A2_d)*Y(:,k-3) + kron(I_N,B)* ( u(:,k)...
                                            + 0.004*[abs(Y(1,k)+1) - abs(Y(1,k)-1); abs(Y(3,k)+1) - abs(Y(3,k)-1); abs(Y(5,k)+1) - abs(Y(5,k)-1);...
                                                     abs(Y(7,k)+1) - abs(Y(7,k)-1); abs(Y(9,k)+1) - abs(Y(9,k)-1) ] ) + ...
                            t*[ 0;  -0.01*cos(Y(1,k));  0;  -0.01*cos(Y(3,k));  0;  -0.01*cos(Y(5,k));  0;  -0.01*cos(Y(7,k)); 0;  -0.01*cos(Y(9,k)) ] );
     elseif k == 3
     Y(:,k+1) =  h1(:,k) * (  kron(I_N,A1)*Y(:,k) + kron(I_N,B)* ( u(:,k)...
                                            + 0.004*[abs(Y(1,k)+1) - abs(Y(1,k)-1); abs(Y(3,k)+1) - abs(Y(3,k)-1); abs(Y(5,k)+1) - abs(Y(5,k)-1); ...
                                                     abs(Y(7,k)+1) - abs(Y(7,k)-1); abs(Y(9,k)+1) - abs(Y(9,k)-1) ] )  + ...
                            t*[ 0;  -0.01*cos(Y(1,k));  0;  -0.01*cos(Y(3,k));  0;  -0.01*cos(Y(5,k));  0;  -0.01*cos(Y(7,k)); 0;  -0.01*cos(Y(9,k)) ] ) +...
                 h2(:,k) * (  kron(I_N,A2)*Y(:,k) + kron(I_N,B)* ( u(:,k)...
                                            + 0.004*[abs(Y(1,k)+1) - abs(Y(1,k)-1); abs(Y(3,k)+1) - abs(Y(3,k)-1); abs(Y(5,k)+1) - abs(Y(5,k)-1);...
                                                     abs(Y(7,k)+1) - abs(Y(7,k)-1); abs(Y(9,k)+1) - abs(Y(9,k)-1) ] ) + ...
                            t*[ 0;  -0.01*cos(Y(1,k));  0;  -0.01*cos(Y(3,k));  0;  -0.01*cos(Y(5,k));  0;  -0.01*cos(Y(7,k)); 0;  -0.01*cos(Y(9,k)) ] );
     elseif k == 2
     Y(:,k+1) =  h1(:,k) * (  kron(I_N,A1)*Y(:,k) + kron(I_N,B)* ( u(:,k)...
                                            + 0.004*[abs(Y(1,k)+1) - abs(Y(1,k)-1); abs(Y(3,k)+1) - abs(Y(3,k)-1); abs(Y(5,k)+1) - abs(Y(5,k)-1); ...
                                                     abs(Y(7,k)+1) - abs(Y(7,k)-1); abs(Y(9,k)+1) - abs(Y(9,k)-1) ] )  + ...
                            t*[ 0;  -0.01*cos(Y(1,k));  0;  -0.01*cos(Y(3,k));  0;  -0.01*cos(Y(5,k));  0;  -0.01*cos(Y(7,k)); 0;  -0.01*cos(Y(9,k)) ] ) +...
                 h2(:,k) * (  kron(I_N,A2)*Y(:,k) + kron(I_N,B)* ( u(:,k)...
                                            + 0.004*[abs(Y(1,k)+1) - abs(Y(1,k)-1); abs(Y(3,k)+1) - abs(Y(3,k)-1); abs(Y(5,k)+1) - abs(Y(5,k)-1);...
                                                     abs(Y(7,k)+1) - abs(Y(7,k)-1); abs(Y(9,k)+1) - abs(Y(9,k)-1) ] ) + ...
                            t*[ 0;  -0.01*cos(Y(1,k));  0;  -0.01*cos(Y(3,k));  0;  -0.01*cos(Y(5,k));  0;  -0.01*cos(Y(7,k)); 0;  -0.01*cos(Y(9,k)) ] );
     else
     Y(:,k+1) =  h1(:,k) * (  kron(I_N,A1)*Y(:,k) + kron(I_N,B)* ( u(:,k)...
                                            + 0.004*[abs(Y(1,k)+1) - abs(Y(1,k)-1); abs(Y(3,k)+1) - abs(Y(3,k)-1); abs(Y(5,k)+1) - abs(Y(5,k)-1); ...
                                                     abs(Y(7,k)+1) - abs(Y(7,k)-1); abs(Y(9,k)+1) - abs(Y(9,k)-1) ] )  + ...
                            t*[ 0;  -0.01*cos(Y(1,k));  0;  -0.01*cos(Y(3,k));  0;  -0.01*cos(Y(5,k));  0;  -0.01*cos(Y(7,k)); 0;  -0.01*cos(Y(9,k)) ] ) +...
                 h2(:,k) * (  kron(I_N,A2)*Y(:,k) + kron(I_N,B)* ( u(:,k)...
                                            + 0.004*[abs(Y(1,k)+1) - abs(Y(1,k)-1); abs(Y(3,k)+1) - abs(Y(3,k)-1); abs(Y(5,k)+1) - abs(Y(5,k)-1);...
                                                     abs(Y(7,k)+1) - abs(Y(7,k)-1); abs(Y(9,k)+1) - abs(Y(9,k)-1) ] ) + ...
                            t*[ 0;  -0.01*cos(Y(1,k));  0;  -0.01*cos(Y(3,k));  0;  -0.01*cos(Y(5,k));  0;  -0.01*cos(Y(7,k)); 0;  -0.01*cos(Y(9,k)) ] );
     end
     
 end
%% 
figure(1);
plot(time,y_R1(1,:)/100,'k-',time,Tr1(1,:)/100,time,Tr1(2,:)/100,time,Tr1(3,:)/100,time,Tr1(4,:)/100,time,Tr1(5,:)/100,'LineWidth',1);
xlabel('time(s)','FontName','Times New Roman','FontSize',12);
ylabel('output response','FontName','Times New Roman','FontSize',12);
f1=legend('$r_{k}$','agent1','agent2','agent3','agent4','agent5','Location','NorthWest','FontName','Times New Roman');
set(f1,'Interpreter','latex','fontsize',10);

figure(3);
plot(time,u(1,:)/100,time,u(2,:)/100,time,u(3,:)/100,time,u(4,:)/100,time,u(5,:)/100,'LineWidth',1);
xlabel('time(s)','FontName','Times New Roman','FontSize',12);
ylabel('control input','FontName','Times New Roman','FontSize',12);
legend('agent1','agent2','agent3','agent4','agent5','FontName','Times New Roman','FontSize',10);

figure(4);
plot(time,e_PID(1,:)/100,time,e_PID(2,:)/100,time,e_PID(3,:)/100,time,e_PID(4,:)/100,time,e_PID(5,:)/100,'LineWidth',1);
xlabel('time(s)','FontName','Times New Roman','FontSize',12);
ylabel('tracking error','FontName','Times New Roman','FontSize',12);
legend('agent1','agent2','agent3','agent4','agent5','FontName','Times New Roman','FontSize',10);