#Dados acoes do dia 25/08/2022

IBOV = [-6.57;-6.74;-1.53;2.85;6.98;0.89;6.06;-10.10;3.22;-11.50;4.69;6.16];
BBAS = [-5.03;-1.38;11.72;-9.39;13.21;7.81;-1.45;-5.94;12.19;-8.85;7.76;15.90];
CPLE = [7.36;-11.67;0.55;8.99;10.85;0.20;6.69;-3.75;8.95;-7.98;3.08;0.23];
CIEL = [-20.30;-4.25;-4.43;14.95;0.90;11.56;21.12;9.21;17.17;-5.14;18.16;26.15];
PETR = [0.82;-1.71;8.93;1.86;14.89;3.12;-3.11;-5.11;-0.54;-8.18;21.02;0.57];
RANI = [-10.09;2.56;10.00;-1.99;-6.96;-5.14;11.66;-5.29;15.37;-10.23;3.45;10.29];
BEEF = [25.00;-6.79;-11.91;23.54;-10.00;13.31;17.11;3.71;8.30;-6.68;-1.28;17.94];
CRFB = [-2.15;-6.31;-7.87;-0.52;9.31;14.28;18.22;-8.66;-5.93;-13.85;11.76;6.01];
RAIL = [-10.31;-4.76;9.88;1.08;-12.05;-0.77;19.68;-11.75;8.00;-9.67;9.83;15.22];
MULT = [-9.81;-2.68;5.83;-4.49;14.85;2.74;10.82;0.69;-2.27;-8.88;8.43;1.26];
HYPE = [-9.64;-12.69;-3.35;4.24;10.05;7.68;15.55;-3.33;3.66;-1.93;11.96;1.13];
SELIC_27 = [0.55;0.61;0.69;0.79;0.89;1.13;0.84;0.72;1.14;1.10;1.01;1.16];
MONTH = 1:12;

M = [IBOV,BBAS,CPLE,CIEL,PETR,RANI,BEEF,CRFB,RAIL,MULT,HYPE];
Retorno_medio(1) = sum(SELIC_27)/12;
Risco(1) = 0;
for i = 2:12
  Retorno_medio(i) = sum(M(:,i-1))/12;
  Risco(i) = sqrt(var(M(:,i-1)));
end

