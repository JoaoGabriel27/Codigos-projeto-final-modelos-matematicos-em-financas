#Comparando o modelo de markowitz e o de indice unico
close all
clear all
pkg load optim

run("Dados_acoes.m")

#Calculando os dados necessarios para o modelo de indice unico

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
  e_1(i) = (1/12)*sum((M(:,i+1) - (alpha(i) + beta(i)*M(:,1))).^2);
end

#Calculando as variancas dos ativos
for i = 1:10
  var_ati1(i) = (beta(i)^2)*var_m + e_1(i);
end

#Calculando a matriz de covariancias
for i = 1:10
  for j = i:10
    Q2(i,j) = beta(i)*beta(j)*var_m;
  end
  Q2(i,i) = var_ati1(i);
end
Q2 = Q2' + Q2 - diag(diag(Q2));

#Calculando os dados necessarios para o modelo de markowitz

Q1 = cov(M(:,2:11));

beq = 1;
Aeq = ones(1,10);
cart_otm_vd_mark = quadprog(Q1,[],[],[],Aeq,beq,[],[]);#Otimizador com venda a descoberto
Retorno_otm_vd_mark = Retorno_medio(3:12)*cart_otm_vd_mark;
Risco_otm_vd_mark = sqrt(cart_otm_vd_mark'*Q1*cart_otm_vd_mark);
Retorno_otm_vd_meses_mark = M(:,2:11)*cart_otm_vd_mark;

cart_otm_vd_ind_uni = quadprog(Q2,[],[],[],Aeq,beq,[],[]);#Otimizador com venda a descoberto
Retorno_otm_vd_ind_uni = Retorno_medio(3:12)*cart_otm_vd_ind_uni;
Risco_otm_vd_ind_uni = sqrt(cart_otm_vd_ind_uni'*Q2*cart_otm_vd_ind_uni);
Retorno_otm_vd_meses_ind_uni = M(:,2:11)*cart_otm_vd_ind_uni;

lb1 = zeros(10,1);
cart_otm_mark = quadprog(Q1,[],[],[],Aeq,beq,lb1,[]);#Otimizador sem venda a descoberto
Retorno_otm_mark = Retorno_medio(3:12)*cart_otm_mark;
Risco_otm_mark = sqrt(cart_otm_mark'*Q1*cart_otm_mark);
Retorno_otm_meses_mark = M(:,2:11)*cart_otm_mark;

cart_otm_ind_uni = quadprog(Q2,[],[],[],Aeq,beq,lb1,[]);#Otimizador sem venda a descoberto
Retorno_otm_ind_uni = Retorno_medio(3:12)*cart_otm_ind_uni;
Risco_otm_ind_uni = sqrt(cart_otm_ind_uni'*Q2*cart_otm_ind_uni);
Retorno_otm_meses_ind_uni = M(:,2:11)*cart_otm_ind_uni;

for i = 1:12
  Retorno_otm_vd_acu_mark(i) = sum(Retorno_otm_vd_meses_mark(1:i));
  Retorno_otm_acu_mark(i) = sum(Retorno_otm_meses_mark(1:i));
  Retorno_otm_vd_acu_ind_uni(i) = sum(Retorno_otm_vd_meses_ind_uni(1:i));
  Retorno_otm_acu_ind_uni(i) = sum(Retorno_otm_meses_ind_uni(1:i));
  Retorno_acu_mercado(i) = sum(M(1:i,1));
end

%Fronteira eficiente com venda a descoberto
A1 = -Retorno_medio(3:12);
R1max = max(Retorno_medio);
Ris1max = max(Risco);
R = 0;
c = 0;
while R < Ris1max+5 #Fronteira eficiente para o caso com venda a descoberto
  b = -(Retorno_otm_vd_mark + c*0.05);
  z1 = quadprog(Q1,[],A1,b,Aeq,beq,[],[]);
  R = sqrt(z1'*Q1*z1);
  X1_mark(c+1,:) = [Retorno_medio(3:12)*z1 R];
  c = c + 1;
end
R = 0;
c = 0;
while R < Ris1max+5 #Fronteira eficiente para o caso com venda a descoberto
  b = -(Retorno_otm_vd_ind_uni + c*0.05);
  z1 = quadprog(Q2,[],A1,b,Aeq,beq,[],[]);
  R = sqrt(z1'*Q2*z1);
  X1_ind_uni(c+1,:) = [Retorno_medio(3:12)*z1 R];
  c = c + 1;
end

%Fronteira eficiente sem venda a descoberto
N = 500;
for i = 1:N#Fronteira eficiente para o caso sem venda a descoberto
  b = -(Retorno_otm_mark + (i-1)*((R1max - Retorno_otm_mark)/(N-1)));
  z2 = quadprog(Q1,[],A1,b,Aeq,beq,lb1,[]);
  X2_mark(i,:) = [Retorno_medio(3:12)*z2 sqrt(z2'*Q1*z2)];
end

for i = 1:N#Fronteira eficiente para o caso sem venda a descoberto
  b = -(Retorno_otm_ind_uni + (i-1)*((R1max - Retorno_otm_ind_uni)/(N-1)));
  z2 = quadprog(Q2,[],A1,b,Aeq,beq,lb1,[]);
  X2_ind_uni(i,:) = [Retorno_medio(3:12)*z2 sqrt(z2'*Q2*z2)];
end

#Grafico comparando os retornos acumulados das carteiras de minimo risco com venda a descoberto dos modelos de indice unico e markowitz
figure(1)
plot(1:12,Retorno_otm_vd_acu_mark,"b","linewidth",2)
hold on
plot(1:12,Retorno_otm_vd_acu_ind_uni,"r","linewidth",2)
plot(1:12,Retorno_acu_mercado,"k","linewidth",2)
xticks([1:1:12])
yticks([-25:5:50])
xlim([0.8 12.2])
ylim([-20 45])
xlabel('Meses')
ylabel('Retorno acumulado (%)')
title("Comparação retorno mensal acumulado das carteira de mínimo risco por Markowitz e por Índ. Úni. com venda a descoberto")
legend({"x_{S}","x_{IU_{S}}","IBOV"},'Location','south','NumColumns',12)
grid on
SizeFontTick = 25;
set(gca,'FontSize',SizeFontTick)

#Grafico comparando os retornos acumulados das carteiras de minimo risco sem venda a descoberto dos modelos de indice unico e markowitz
figure(2)
plot(1:12,Retorno_otm_acu_mark,"b","linewidth",2)
hold on
plot(1:12,Retorno_otm_acu_ind_uni,"r","linewidth",2)
plot(1:12,Retorno_acu_mercado,"k","linewidth",2)
xticks([1:1:12])
yticks([-25:5:50])
xlim([0.8 12.2])
ylim([-20 40])
xlabel('Meses')
ylabel('Retorno acumulado (%)')
title("Comparação retorno mensal acumulado das carteira de mínimo risco por Markowitz e por Índ. Úni. com venda a descoberto")
legend({"x_{NS}","x_{IU_{NS}}","IBOV"},'Location','south','NumColumns',12)
grid on
SizeFontTick = 25;
set(gca,'FontSize',SizeFontTick)

#Grafico comparando as fronteiras eficientes com venda a descoberto dos modelos de indice unico e markowitz
figure(3)
plot(X1_mark(:,2),X1_mark(:,1),"b","linewidth",2)
hold on
plot(X1_ind_uni(:,2),X1_ind_uni(:,1),"r","linewidth",2)
scatter(Risco_otm_vd_mark,Retorno_otm_vd_mark,'filled')
scatter(Risco_otm_vd_ind_uni,Retorno_otm_vd_ind_uni,'filled')
scatter(Risco(2),Retorno_medio(2),'filled');
str1 = {"IBOV"};
xticks([0:0.5:50])
yticks([-10:1:50])
xlim([1.3 14.1])
ylim([-0.9 18.1])
xlabel('Risco (%)')
ylabel('Retorno médio (%)')
title("Comparação carteiras e fronteiras eficientes no plano Risco x Retorno, permitida venda a descoberto")
legend("Fronteira Eficiente por Markowitz","Fronteira Eficiente por índice único","x_{S}","x_{IU_{S}}",'Location','northwest')
text(Risco(2).-0.2,Retorno_medio(2).+0.4,"IBOV",'FontSize',20)
grid on
SizeFontTick = 25;
set(gca,'FontSize',SizeFontTick)

#Grafico comparando as fronteiras eficientes sem venda a descoberto dos modelos de indice unico e markowitz
figure(4)
plot(X2_mark(:,2),X2_mark(:,1),"b","linewidth",2)
hold on
plot(X2_ind_uni(:,2),X2_ind_uni(:,1),"r","linewidth",2)
scatter(Risco_otm_mark,Retorno_otm_mark,'filled')
scatter(Risco_otm_ind_uni,Retorno_otm_ind_uni,'filled')
scatter(Risco(2),Retorno_medio(2),'filled');
str1 = {"IBOV"};
xticks([0:0.5:50])
yticks([-10:0.5:50])
xlim([4.5 13.9])
ylim([-0.9 8.1])
xlabel('Risco (%)')
ylabel('Retorno médio (%)')
title("Comparação carteiras e fronteiras eficientes no plano Risco x Retorno, proibida venda a descoberto")
legend("Fronteira Eficiente por Markowitz","Fronteira Eficiente por índice único","x_{NS}","x_{IU_{NS}}",'Location','northwest')
text(Risco(2).-0.2,Retorno_medio(2).+0.4,"IBOV",'FontSize',20)
grid on
SizeFontTick = 25;
set(gca,'FontSize',SizeFontTick)


