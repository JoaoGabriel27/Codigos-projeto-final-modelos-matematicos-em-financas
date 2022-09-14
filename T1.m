%
#Parte T1 Projeto final modelos de matematica em financas
close all
clear all
pkg load optim

run("Dados_acoes.m")

#Parte (a)
%}
Q = cov(M(:,2:11));#Matriz de covariancia
%
beq = 1;
Aeq = ones(1,10);
cart_otm_vd = quadprog(Q,[],[],[],Aeq,beq,[],[]);#Otimizador com venda a descoberto
Retorno_otm_vd = Retorno_medio(3:12)*cart_otm_vd;#Retorno medio da carteira otima com venda a descoberto
Risco_otm_vd = sqrt(cart_otm_vd'*Q*cart_otm_vd);#Risco da carteira otima com venda a descoberto
Retorno_otm_vd_meses = M(:,2:11)*cart_otm_vd;#Retorno mes-a-mes da carteira otima com venda a descoberto
%{
#Grafico flutuação do retorno
figure(1)
plot(MONTH,Retorno_otm_vd_meses,"k","linewidth",2)
hold on
cores = ['#FF0000';'#00FF00';'#0000FF';'#00FFFF';'#FF00FF';'#FFFF00';'#D95319';'#EDB120';'#7E2F8E';'#77AC30';'#A2142F'];
for i = 1:11
  plot(MONTH,M(:,i),"linewidth",1,"color",cores(i,:))
end
xticks([1:1:12])
yticks([-40:2:40])
xlim([0.8 12.2])
ylim([-22.5 34.5])
xlabel('Meses')
ylabel('Retorno (%)')
title("Variação mensal do retorno dos ativos e da carteira de mínimo risco com venda a descoberto")
legend({"Otm vd","IBOV","BBAS","CPLE","CIEL","PETR","RANI","BEEF","CRFB","RAIL","MULT","HYPE"},'Location','south','NumColumns',12)
grid on
SizeFontTick = 20;
set(gca,'FontSize',SizeFontTick)
%}
#Grafico retorno acumulado
for i = 1:12
  Retorno_otm_vd_acu(i) = sum(Retorno_otm_vd_meses(1:i));#Retornos acumulados de cada mes da carteira otima com venda a descoberto
end
M_acu(1,:) = M(1,:);#Retornos acumulados das acoes
for i = 2:12
  M_acu(i,:) = sum(M(1:i,:));
end
%
figure(2)
plot(MONTH,Retorno_otm_vd_acu,"k","linewidth",2)
hold on
cores = ['#FF0000';'#00FF00';'#0000FF';'#00FFFF';'#FF00FF';'#FFFF00';'#D95319';'#EDB120';'#7E2F8E';'#77AC30';'#A2142F'];
for i = 1:11
  plot(MONTH,M_acu(:,i),"linewidth",1,"color",cores(i,:))
end
xticks([1:1:12])
yticks([-40:5:100])
xlim([0.8 12.2])
ylim([-40 90])
xlabel('Meses')
ylabel('Retorno acumulado (%)')
title("Retorno mensal acumulado dos ativos e da carteira de mínimo risco com venda a descoberto")
legend({"x_{S}","IBOV","BBAS3","CPLE11","CIEL3","PETR3","RANI3","BEEF3","CRFB3","RAIL3","MULT3","HYPE3"},'Location','south','NumColumns',12)
grid on
SizeFontTick = 25;
set(gca,'FontSize',SizeFontTick)
%

%Fronteira eficiente com venda a descoberto
A1 = -Retorno_medio(3:12);
R1max = max(Retorno_medio);
Ris1max = max(Risco);
R = 0;
c = 0;
while R < Ris1max+5 #Fronteira eficiente para o caso com venda a descoberto
  b = -(Retorno_otm_vd + c*0.05);
  z1 = quadprog(Q,[],A1,b,Aeq,beq,[],[]);
  R = sqrt(z1'*Q*z1);
  X1(c+1,:) = [Retorno_medio(3:12)*z1 R];
  c = c + 1;
end
%
figure(3)
plot(X1(:,2),X1(:,1),"k","linewidth",2)
hold on
scatter(Risco_otm_vd,Retorno_otm_vd,'filled')
scatter(Risco(2:12),Retorno_medio(2:12),'filled');
str1 = {"IBOV","BBAS3","CPLE11","CIEL3","PETR3","RANI3","BEEF3","CRFB3","RAIL3","MULT3","HYPE3"};
xticks([0:0.5:50])
yticks([-10:1:50])
xlim([1.3 14.1])
ylim([-0.9 18.1])
xlabel('Risco (%)')
ylabel('Retorno médio (%)')
title("Carteiras e fronteira eficiente no plano Risco x Retorno, permitida venda a descoberto")
legend("Fronteira Eficiente","Carteira Ótima","100% no Ativo",'Location','northwest')
text(Risco(2:12).-0.2,Retorno_medio(2:12).+0.4,str1,'FontSize',20)
grid on
SizeFontTick = 25;
set(gca,'FontSize',SizeFontTick)
%}
#Parte (b)

lb1 = zeros(10,1);
cart_otm = quadprog(Q,[],[],[],Aeq,beq,lb1,[]);#Otimizador sem venda a descoberto
Retorno_otm = Retorno_medio(3:12)*cart_otm;#Retorno medio da carteira otima sem venda a descoberto
Risco_otm = sqrt(cart_otm'*Q*cart_otm);#Risco da carteira otima sem venda a descoberto
Retorno_otm_meses = M(:,2:11)*cart_otm;#Retorno mes-a-mes da carteira otima sem venda a descoberto
%{
figure(4)
plot(MONTH,Retorno_otm_meses,"k","linewidth",2)
hold on
for i = 1:11
  plot(MONTH,M(:,i),"linewidth",1,"color",cores(i,:))
end
legend("Otm","IBOV","BBAS","CPLE","CIEL","PETR","RANI","BEEF","CRFB","RAIL","MULT","HYPE",'Location','south','NumColumns',12)
xticks([1:1:12])
yticks([-40:2:40])
xlim([0.8 12.2])
ylim([-22.5 34.5])
xlabel('Meses')
ylabel('Retorno (%)')
title("Variação mensal do retorno dos ativos e da carteira de mínimo risco sem venda a descoberto")
grid on
SizeFontTick = 20;
set(gca,'FontSize',SizeFontTick)
%}
#Grafico retorno acumulado
for i = 1:12
  Retorno_otm_acu(i) = sum(Retorno_otm_meses(1:i));#Retornos acumulados de cada mes da carteira otima sem venda a descoberto
end
%
#Grafico retorno acumulado sem venda a descoberto
figure(5)
plot(MONTH,Retorno_otm_acu,"k","linewidth",2)
hold on
cores = ['#FF0000';'#00FF00';'#0000FF';'#00FFFF';'#FF00FF';'#FFFF00';'#D95319';'#EDB120';'#7E2F8E';'#77AC30';'#A2142F'];
for i = 1:11
  plot(MONTH,M_acu(:,i),"linewidth",1,"color",cores(i,:))
end
xticks([1:1:12])
yticks([-40:5:100])
xlim([0.8 12.2])
ylim([-40 90])
xlabel('Meses')
ylabel('Retorno acumulado (%)')
title("Retorno mensal acumulado dos ativos e da carteira de mínimo risco sem venda a descoberto")
legend({"x_{NS}","IBOV","BBAS3","CPLE11","CIEL3","PETR3","RANI3","BEEF3","CRFB3","RAIL3","MULT3","HYPE3"},'Location','south','NumColumns',12)
grid on
SizeFontTick = 25;
set(gca,'FontSize',SizeFontTick)
%

%Fronteira eficiente sem venda a descoberto
N = 500;
for i = 1:N#Fronteira eficiente para o caso sem venda a descoberto
  b = -(Retorno_otm + (i-1)*((R1max - Retorno_otm)/(N-1)));
  z2 = quadprog(Q,[],A1,b,Aeq,beq,lb1,[]);
  X2(i,:) = [Retorno_medio(3:12)*z2 sqrt(z2'*Q*z2)];
end
%
figure(6)
plot(X2(:,2),X2(:,1),"k","linewidth",2)
hold on
scatter(Risco_otm,Retorno_otm,'filled')
scatter(Risco(2:12),Retorno_medio(2:12),'filled');
str1 = {"IBOV","BBAS3","CPLE11","CIEL3","PETR3","RANI3","BEEF3","CRFB3","RAIL3","MULT3","HYPE3"};
xticks([0:0.5:50])
yticks([-10:0.5:50])
xlim([4.8 14.1])
ylim([-0.6 7.6])
xlabel('Risco (%)')
ylabel('Retorno médio (%)')
title("Carteiras e fronteira eficiente no plano Risco x Retorno, não permitida venda a descoberto")
legend("Fronteira Eficiente","Carteira Ótima","100% no Ativo",'Location','northwest')
text(Risco(2:12).-0.15,Retorno_medio(2:12).+0.2,str1,'FontSize',20)
grid on
SizeFontTick = 25;
set(gca,'FontSize',SizeFontTick)

