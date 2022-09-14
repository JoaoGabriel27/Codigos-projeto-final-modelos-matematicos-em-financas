%
#Parte T2 Projeto final modelos de matematica em financas
close all
clear all
pkg load optim

run("Dados_acoes.m")

R_a = Retorno_medio(3:12);
R_f = Retorno_medio(1);
Q = cov(M(:,2:11));
f = @(x)(-(R_a*x - R_f)/(sqrt(x'*Q*x)));

#Resolvendo max teta(x) sujeito a restriçoes sem venda a descoberto
lb = zeros(10,1);
beq = 1;
Aeq = ones(1,10);
x_0 = ones(10,1)/10;
x_svd = fmincon(f,x_0,[],[],Aeq,beq,lb,[]);#Carteira tangente sem venda a descoberto
teta_svd = -f(x_svd);#Angulo otimo sem venda a descoberto
Reta_svd = @(x)(teta_svd*x + R_f);#Reta tangente sem venda a descoberto

cart_otm_svd = quadprog(Q,[],[],[],Aeq,beq,lb,[]);#Otimizador sem venda a descoberto
Retorno_otm_svd = Retorno_medio(3:12)*cart_otm_svd;#Retorno medio da carteira otima sem venda a descoberto
%Fronteira eficiente com venda a descoberto
A1 = -Retorno_medio(3:12);
R1max = max(Retorno_medio);
N = 500;
for i = 1:N#Fronteira eficiente para o caso sem venda a descoberto
  b = -(Retorno_otm_svd + (i-1)*((R1max - Retorno_otm_svd)/(N-1)));
  z2 = quadprog(Q,[],A1,b,Aeq,beq,lb,[]);
  X2(i,:) = [Retorno_medio(3:12)*z2 sqrt(z2'*Q*z2)];
end

risco_svd = sqrt(x_svd'*Q*x_svd);#Risco da carteira tangente sem venda a descoberto
retorno_svd = Retorno_medio(3:12)*x_svd;#Retorno medio da carteira tangente sem venda a descoberto
x_svd_reta = [0 sqrt(x_svd'*Q*x_svd)];
#Grafico reta tangente sem venda a descoberto
figure(1)
plot(X2(:,2),X2(:,1),"r","linewidth",2)
hold on
plot(x_svd_reta,Reta_svd(x_svd_reta),"b","linewidth",2)
scatter(x_svd_reta(2),Reta_svd(x_svd_reta(2)),'filled')
scatter(0,Retorno_medio(1),'filled');
xticks([0:0.5:50])
yticks([-10:1:50])
xlabel('Risco (%)')
ylabel('Retorno médio (%)')
title("Reta tangente e fronteira eficiente sem venda a descoberto")
legend("Fronteira eficiente","Reta tangente","Carteira tangente","Ativo livre de risco",'Location','northwest')
grid on
SizeFontTick = 25;
set(gca,'FontSize',SizeFontTick)

#Resolvendo max teta(x) sujeito a restriçoes com venda a descoberto
z = Q\(R_a' - R_f);#Por sistema linear
x_cvd = z/sum(z);#Carteira tangente com venda a descoberto
teta_cvd = -f(x_cvd);#Angulo otimo com venda a descoberto
Reta_cvd = @(x)(teta_cvd*x + R_f)#Reta tangente com venda a descoberto

cart_otm_vd = quadprog(Q,[],[],[],Aeq,beq,[],[]);#Otimizador com venda a descoberto
Retorno_otm_vd = Retorno_medio(3:12)*cart_otm_vd;#Retorno medio da carteira otima com venda a descoberto
%Fronteira eficiente com venda a descoberto
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

risco_cvd = sqrt(x_cvd'*Q*x_cvd);#Risco da carteira tangente com venda a descoberto
retorno_cvd = Retorno_medio(3:12)*x_cvd;#Retorno medio da carteira tangente com venda a descoberto
x_cvd_reta = [0 sqrt(x_cvd'*Q*x_cvd) 15];
#Grafico reta tangente com venda a descoberto
figure(2)
plot(X1(:,2),X1(:,1),"r","linewidth",2)
hold on
plot(x_cvd_reta,Reta_cvd(x_cvd_reta),"b","linewidth",2)
scatter(x_cvd_reta(2),Reta_cvd(x_cvd_reta(2)),'filled')
scatter(0,Retorno_medio(1),'filled');
xticks([0:0.5:50])
yticks([-10:1:50])
xlabel('Risco (%)')
ylabel('Retorno médio (%)')
title("Reta tangente e fronteira eficiente com venda a descoberto")
legend("Fronteira eficiente","Reta tangente","Carteira tangente","Ativo livre de risco",'Location','northwest')
grid on
SizeFontTick = 25;
set(gca,'FontSize',SizeFontTick)

#Grafico retorno acumulado
Retorno_svd_meses = M(:,2:11)*x_svd;#Retorno mes-a-mes da carteira tangente sem venda a descoberto
for i = 1:12
  Retorno_acu_svd(i) = sum(Retorno_svd_meses(1:i));#Retorno acumulado da carteira tangente sem venda a descoberto
end
M_acu(1,:) = M(1,:);
for i = 2:12
  M_acu(i,:) = sum(M(1:i,:));
end
%
#Grafico retorno acumulado carteria tangente sem venda a descoberto
figure(3)
plot(1:12,Retorno_acu_svd,"k","linewidth",2)
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
title("Retorno mensal acumulado dos ativos e da carteira tangente sem venda a descoberto")
legend({"T_{NS}","IBOV","BBAS3","CPLE11","CIEL3","PETR3","RANI3","BEEF3","CRFB3","RAIL3","MULT3","HYPE3"},'Location','south','NumColumns',12)
grid on
SizeFontTick = 25;
set(gca,'FontSize',SizeFontTick)

Retorno_cvd_meses = M(:,2:11)*x_cvd;#Retorno mes-a-mes da carteira tangente com venda a descoberto
for i = 1:12
  Retorno_acu_cvd(i) = sum(Retorno_cvd_meses(1:i));#Retorno acumulado da carteira tangente com venda a descoberto
end
%
#Grafico retorno acumulado carteria tangente com venda a descoberto
figure(4)
plot(1:12,Retorno_acu_cvd,"k","linewidth",2)
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
title("Retorno mensal acumulado dos ativos e da carteira tangente com venda a descoberto")
legend({"T_{S}","IBOV","BBAS3","CPLE11","CIEL3","PETR3","RANI3","BEEF3","CRFB3","RAIL3","MULT3","HYPE3"},'Location','south','NumColumns',12)
grid on
SizeFontTick = 25;
set(gca,'FontSize',SizeFontTick)
