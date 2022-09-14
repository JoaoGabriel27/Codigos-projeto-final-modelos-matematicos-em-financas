#figura fronteira eficiente

close all
clear all
pkg load optim
run("Dados_acoes.m")

Q = cov(M(:,2:11));

risco = @(x)(sqrt(x'*Q*x));#Funcao do risco
retorno = @(x)(Retorno_medio(3:12)*x);#Funcao do retorno

#Criando as carteiras "aleatorias"
for i = 1:100
  x(:,i) = rand(10,1);
  x(:,i) = x(:,i)/sum(x(:,i));
end
for i = 100:10:1100
  for j = 1:10
    for k = 1:0.1:2
      x(:,i+j) = rand(10,1);
      x(j,i+j) = j*0.3*k;
      x(:,i+j) = x(:,i+j)/sum(x(:,i+j));
    end
  end
end

#Grafico carteiras aleatorias no plano risco x retorno
figure(1)
scatter(risco(x(:,1)),retorno(x(:,1)),'filled')
hold on
for i = 2:1100
  scatter(risco(x(:,i)),retorno(x(:,i)),'filled')
end
xticks([5:0.2:10])
yticks([1:0.2:5.8])
xlim([5 10])
ylim([1 5.8])
xlabel('Risco (%)')
ylabel('Retorno médio (%)')
title("Carteiras aleatórias")
grid on
SizeFontTick = 20;
set(gca,'FontSize',SizeFontTick)

#Para criar a fronteira eficiente
beq = 1;
Aeq = ones(1,10);
lb1 = zeros(10,1);
cart_otm = quadprog(Q,[],[],[],Aeq,beq,lb1,[]);#Otimizador sem venda a descoberto
Retorno_otm = Retorno_medio(3:12)*cart_otm;
Risco_otm = sqrt(cart_otm'*Q*cart_otm);
Retorno_otm_meses = M(:,2:11)*cart_otm;
%Fronteira eficiente sem venda a descoberto
A1 = -Retorno_medio(3:12);
R1max = max(Retorno_medio);
N = 500;
for i = 1:N#Fronteira eficiente para o caso sem venda a descoberto
  b = -(Retorno_otm + (i-1)*((R1max - Retorno_otm)/(N-1)));
  z2 = quadprog(Q,[],A1,b,Aeq,beq,lb1,[]);
  X2(i,:) = [Retorno_medio(3:12)*z2 sqrt(z2'*Q*z2)];
end

#Grafico das carteiras aleatorias junto da fronteira eficiente no plano risco x retorno
figure(2)
plot(X2(:,2),X2(:,1),"k","linewidth",2)
hold on
for i = 1:1100
  scatter(risco(x(:,i)),retorno(x(:,i)),'filled')
end
xticks([4.6:0.2:10])
yticks([1:0.2:6.6])
xlim([4.6 10])
ylim([1 6.6])
xlabel('Risco (%)')
ylabel('Retorno médio (%)')
title("Fronteira eficiente e carteiras aleatórias")
grid on
SizeFontTick = 20;
set(gca,'FontSize',SizeFontTick)
