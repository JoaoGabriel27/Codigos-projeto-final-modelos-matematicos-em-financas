%
#Parte T3 Projeto final modelos de matematica em financas
close all
clear all
pkg load optim

run("Dados_acoes.m")

#Calculando os betas
var_m = var(M(:,1));
for i = 1:10
  covariancias(i) = cov(M(:,1),M(:,i+1));
  beta(i) = cov(M(:,1),M(:,i+1))/var_m;
end

#Calculando os retornos medios
Retorno_medio_new = Retorno_medio(2:12);

#Calculando os alphas
for i = 1:10
  alpha(i) = Retorno_medio_new(i+1) - beta(i)*Retorno_medio_new(1);
end

#Calculando as variancas do erro
for i = 1:10
  e_1(i) = (1/12)*sum((M(:,i+1) - (alpha(i) + beta(i)*M(:,1))).^2);#Valor usado pois assim foi pedido
  e_2(i) = var(M(:,i+1)) - (beta(i)^2)*var_m;#Outro modo de calcular
  e_3(:,i) = M(:,i+1) - alpha(i) - beta(i)*M(:,1);#Outro modo de calcular
end

#Calculando as variancas dos ativos
for i = 1:10
  var_ati1(i) = (beta(i)^2)*var_m + e_1(i);
end

#Calculando a matriz de covariancias
for i = 1:10
  for j = i:10
    Q(i,j) = beta(i)*beta(j)*var_m;
  end
  Q(i,i) = var_ati1(i);
end
Q = Q' + Q - diag(diag(Q));

#Calculando a carteira de minimo risco

#Com venda a descoberto
beq = 1;
Aeq = ones(1,10);
cart_otm_vd = quadprog(Q,[],[],[],Aeq,beq,[],[]);#Otimizador com venda a descoberto
Retorno_otm_vd = Retorno_medio_new(2:11)*cart_otm_vd;#Retorno medio carteira otima com venda a descoberto
Risco_otm_vd = sqrt(cart_otm_vd'*Q*cart_otm_vd);#Risco carteira otima com venda a descoberto
Retorno_otm_vd_meses = M(:,2:11)*cart_otm_vd;#Retorno mes-a-mes carteira otima com venda a descoberto

#Sem venda a descoberto
lb1 = zeros(10,1);
cart_otm = quadprog(Q,[],[],[],Aeq,beq,lb1,[]);#Otimizador sem venda a descoberto
Retorno_otm = Retorno_medio_new(2:11)*cart_otm;#Retorno medio carteira otima sem venda a descoberto
Risco_otm = sqrt(cart_otm'*Q*cart_otm);#Risco carteira otima sem venda a descoberto
Retorno_otm_meses = M(:,2:11)*cart_otm;#Retorno mes-a-mes carteira otima sem venda a descoberto

for i = 1:12
  Retorno_otm_vd_acu(i) = sum(Retorno_otm_vd_meses(1:i));#Retorno acumulado carteira otima com venda a descoberto
end
M_acu(1,:) = M(1,:);#Retorno acumulado ativos
for i = 2:12
  M_acu(i,:) = sum(M(1:i,:));
end
%
#Grafico retorno acumulado carteria otima com venda a descoberto
figure(1)
plot(1:12,Retorno_otm_vd_acu,"k","linewidth",2)
hold on
cores = ['#FF0000';'#00FF00';'#0000FF';'#00FFFF';'#FF00FF';'#FFFF00';'#D95319';'#EDB120';'#7E2F8E';'#77AC30';'#A2142F'];
for i = 1:11
  plot(1:12,M_acu(:,i),"linewidth",1,"color",cores(i,:))
end
xticks([1:1:12])
yticks([-40:5:100])
xlim([0.8 12.2])
ylim([-40 90])
xlabel('Meses')
ylabel('Retorno acumulado (%)')
title("Retorno mensal acumulado dos ativos e da carteira de mínimo risco do modelo de índice único com venda a descoberto")
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
#Grafico fronteira eficiente carteira otima com venda a descoberto
figure(2)
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
title("Carteiras e fronteira eficiente do modelo de índice único no plano Risco x Retorno, permitida venda a descoberto")
legend("Fronteira Eficiente","Carteira Ótima","100% no Ativo",'Location','northwest')
text(Risco(2:12).-0.2,Retorno_medio(2:12).+0.4,str1,'FontSize',20)
grid on
SizeFontTick = 25;
set(gca,'FontSize',SizeFontTick)

for i = 1:12
  Retorno_otm_acu(i) = sum(Retorno_otm_meses(1:i));#Retorno acumulado carteria otima sem venda a descoberto
end
%
#Grafico retorno acumulado carteria otima sem venda a descoberto
figure(3)
plot(1:12,Retorno_otm_acu,"k","linewidth",2)
hold on
cores = ['#FF0000';'#00FF00';'#0000FF';'#00FFFF';'#FF00FF';'#FFFF00';'#D95319';'#EDB120';'#7E2F8E';'#77AC30';'#A2142F'];
for i = 1:11
  plot(1:12,M_acu(:,i),"linewidth",1,"color",cores(i,:))
end
xticks([1:1:12])
yticks([-40:5:100])
xlim([0.8 12.2])
ylim([-40 90])
xlabel('Meses')
ylabel('Retorno acumulado (%)')
title("Retorno mensal acumulado dos ativos e da carteira de mínimo risco do modelo de índice único sem venda a descoberto")
legend({"x_{NS}","IBOV","BBAS3","CPLE11","CIEL3","PETR3","RANI3","BEEF3","CRFB3","RAIL3","MULT3","HYPE3"},'Location','south','NumColumns',12)
grid on
SizeFontTick = 25;
set(gca,'FontSize',SizeFontTick)
%

%Fronteira eficiente sem venda a descoberto
N = 500;
for i = 1:N#Fronteira eficiente para o caso sem venda a descoberto
  b = -(Retorno_otm + (i)*((R1max - Retorno_otm)/(N)));
  z2 = quadprog(Q,[],A1,b,Aeq,beq,lb1,[]);
  X2(i,:) = [Retorno_medio(3:12)*z2 sqrt(z2'*Q*z2)];
end
%
#Grafico fronteira eficiente carteira otima sem venda a descoberto
figure(4)
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
title("Carteiras e fronteira eficiente do modelo de índice único no plano Risco x Retorno, não permitida venda a descoberto")
legend("Fronteira Eficiente","Carteira Ótima","100% no Ativo",'Location','northwest')
text(Risco(2:12).-0.15,Retorno_medio(2:12).+0.2,str1,'FontSize',20)
grid on
SizeFontTick = 25;
set(gca,'FontSize',SizeFontTick)



