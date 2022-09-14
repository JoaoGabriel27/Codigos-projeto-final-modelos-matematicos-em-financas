%
#Parte T4 Projeto final modelos de matematica em financas
close all
clear all
pkg load optim

run("Dados_acoes.m")

#Reta do mercado

r = @(x)(Retorno_medio(1) + ((Retorno_medio(2) - Retorno_medio(1))/Risco(2)*x))#Onde x é o risco de uma carteira eficiente

x = [0 Risco(2) Risco(2)+15]

#Reta do mercado e os ativos no plano risco x retorno
figure(1)
plot(x,r(x),"k","linewidth",2)
hold on
for i = 1:2
  scatter(Risco(i),Retorno_medio(i),"filled")
end
scatter(Risco(3:12),Retorno_medio(3:12),"filled")
str1 = {"BBAS3","CPLE11","CIEL3","PETR3","RANI3","BEEF3","CRFB3","RAIL3","MULT3","HYPE3"};
xticks([0:0.5:50])
yticks([-10:1:50])
xlim([0 16.1])
ylim([-4 8.1])
xlabel('Risco (%)')
ylabel('Retorno médio (%)')
title("Reta do mercado e ativos no plano Risco x Retorno")
legend("Reta do Mercado","100% no Ativo Livre de Risco","100% no Mercado","100% nos Ativos com Risco",'Location','northwest')
text(Risco(3:12).-0.2,Retorno_medio(3:12).+0.4,str1,'FontSize',20)
grid on
SizeFontTick = 25;
set(gca,'FontSize',SizeFontTick)

rt = @(x)(Retorno_medio(1) + (Retorno_medio(2) - Retorno_medio(1))*x)#Reta de titulos, onde x é o beta

#Calculando os betas
var_m = var(M(:,1));
beta = [0 1];
for i = 3:12
  beta(i) = cov(M(:,1),M(:,i-1))/var_m;
end

#Reta de titulos com os betas aplicados a ela
figure(2)
plot(x,rt(x),"k","linewidth",2)
hold on
cores = ['#FF0000';'#00FF00';'#0000FF';'#00FFFF';'#FF00FF';'#FFFF00';'#D95319';'#EDB120';'#7E2F8E';'#77AC30';'#A2142F';'#0072BD'];
for i = 1:12
  scatter(beta(i),rt(beta(i)),"MarkerEdgeColor",cores(i,:),"MarkerFaceColor",cores(i,:))
end
str1 = {"Selic27","IBOV","BBAS3","CPLE11","CIEL3","PETR3","RANI3","BEEF3","CRFB3","RAIL3","MULT3","HYPE3"};
xticks([0:0.1:5])
yticks([-10:0.25:10])
xlim([0 1.6])
ylim([-1.6 1.8])
xlabel('\beta',"interpreter","tex")
ylabel('Retorno médio (%)')
title("Reta do mercado aplicando o \beta dos ativos no plano Beta x Retorno",'interpreter', 'tex')
legend("Reta de Títulos","Selic27","IBOV","BBAS3","CPLE11","CIEL3","PETR3","RANI3","BEEF3","CRFB3","RAIL3","MULT3","HYPE3")
grid on
SizeFontTick = 25;
set(gca,'FontSize',SizeFontTick)




