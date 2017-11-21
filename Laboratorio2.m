%% Lab2
%
clear all;
close all;

%% Tarea 1
%
B1=[1 2 1 2 1];
A1=1;
B2=[1 0 0 0 0 -1];
A2=1;
B3=[1 6 11 0 6];
A3=1;
B4=[1 -2 3];
A4=[3 -2 1];
B5=conv([1 -sqrt(2) 1],[1 1]);
A5=conv([1 0.81],[1 0.8*sqrt(2) 0.64]);

tarea1(B1,A1,'\frac{1+2z^{-1}+z^{-2}+2z^{-3}+z^{-4}}{1}');
tarea1(B2,A2,'');
tarea1(B3,A3,'');
tarea1(B4,A4,'');
tarea1(B5,A5,'');

% Retardo(--> orden=mitad retardo)

% NUM: fir/iir (ver polos) -- Resp.Frec -- Fase lineal -- 
% 1  : FIR -- paso bajo -- lineal 
% 2  : FIR -- notch (comb porque son equispaciados) -- lineal
% 3  : FIR -- p bajo -- no lilneal (no hay simetria de coeficientes)
% 4  : IIR -- p todo -- no
% 5  : IIR -- respuesta arbitraria -- no

%% Tarea 2
% per i filtri fir 'H(z)=sum b_k z^{-k}' e i coeff di H(z) corrispondono 
% alla risposta impulsional h(n). Quindi si fa:
% h=fir1(N,wc,ventana)
% B={h(n)}
%
% domanda slide: La atenuacion depende solo de la ventana
clear
close all;
Fm=20e3;
F=5e3;
N=1000; % incrementamos N 20, 1000: se estrecha la banda de transicion

figure
h_rect=fir1(N,F/(Fm/2),rectwin(N+1));
tarea1(h_rect,1,'');

figure
h_hamm=fir1(N,F/(Fm/2),hamming(N+1));
tarea1(h_hamm,1,'');

fvtool(h_rect,1,h_hamm,1);

%% Tarea 3
%
clear
close all
Fm = 1000; % Frecuencia de muestreo
RP = 3; % Rizado en la banda pasante en DB
RS = 60; % Rizado en la banda atenuada en DB sin signo
F1 = 200; % Inicio de la banda de transición
W = [1 0]; % Valores ideales del filtro
r1 = (10^(RP/20)-1)/(10^(RP/20)+1); % Rizado en la banda pasante en escala lineal
r2 = 10^(-RS/20); % Rizado en la banda atenuada en escala lineal
rizado = [r1 r2];
F2 = 300; % Fin de la banda de transición
[n,fo,mo,w] = firpmord( [F1 F2], W, rizado, Fm );
b = firpm(n,fo,mo,w);
% Representación en dB
% FREQZ(B,A,N,Fs),
freqz(b,1,Fm,Fm) % Respresentación con N=Fm (resolución de 1 Hz)
% Representación en escala lineal
[H,w] = freqz(b,1,Fm,Fm);
figure
plot(w,abs(H))
xlabel('Frecuencia (Hz)')
ylabel('|H(\omega)|')
grid on

% looop over F2
Fm = 1000; % Frecuencia de muestreo
RP = 3; % Rizado en la banda pasante en DB
RS = 60; % Rizado en la banda atenuada en DB sin signo
F1 = 200; % Inicio de la banda de transición
F2 = 210:10:300; % Anchura banda transicion
W = [1 0]; % Valores ideales del filtro
r1 = (10^(RP/20)-1)/(10^(RP/20)+1); % Rizado en la banda pasante en escala lineal
r2 = 10^(-RS/20); % Rizado en la banda atenuada en escala lineal
rizado = [r1 r2];
orden = zeros(1,length(F2));
for i=1:length(F2)
    [n,fo,mo,w] = firpmord( [F1 F2(i)], W, rizado, Fm );
    orden(i) = n;
end
figure
plot(F2,orden);
xlabel('Banda transicion (Hz)')
ylabel('Orden')
% el incremento de la atenuacion es lineal con el aumento del orden

%loop over RS
Fm = 1000; % Frecuencia de muestreo
RP = 3; % Rizado en la banda pasante en DB
RS = 20:10:100; % Rizado en la banda atenuada en DB sin signo
F1 = 200; % Inicio de la banda de transición
F2 = 300; % Anchura banda transicion
W = [1 0]; % Valores ideales del filtro
r1 = (10^(RP/20)-1)/(10^(RP/20)+1); % Rizado en la banda pasante en escala lineal
r2 = 10.^(-RS./20); % Rizado en la banda atenuada en escala lineal

orden = zeros(1,length(RS));
for i=1:length(RS)
    rizado = [r1 r2(i)];
    [n,fo,mo,w] = firpmord( [F1 F2], W, rizado, Fm );
    orden(i) = n;
end
figure
plot(RS,orden);
xlabel('Atenuacion')
ylabel('Orden')

%% Tarea 4
% un'onda quadra de freq 100Hz, lo spettro ha solo le frequenze dispari
% multiplo di 100. La filtro con un filtro di Fpassante=200Hz, Fstop=300.
% Lascia passare l'armonica di 100Hz e le altre le attenua. Se
% l'attenuazione è molto alta, l'uscita sarà una sinusoide.
Fm=1000; % Frecuencia de muestreo (Hz)
Fp=200; % Frecuencia limite de la banda de paso (Hz)
Rp=2; % Rizado banda pasante (dB)
Rs=40; % Rizado banda de stop (dB)
BT=100; % Banda de transicion (Hz)

[N_butt,Wn_butt]=buttord(Fp/(Fm/2),(Fp+BT)/(Fm/2),Rp,Rs);
[B_butt,A_butt]=butter(N_butt,Wn_butt);
tarea1(B_butt,A_butt,'Butter filter');

[N_cheb1,Wn_cheb1]=cheb1ord(Fp/(Fm/2),(Fp+BT)/(Fm/2),Rp,Rs);
[B_cheb1,A_cheb1]=cheby1(N_cheb1,Rp,Wn_cheb1);
tarea1(B_cheb1,A_cheb1,'Tchebyshev 1 filter');

[N_cheb2,Wn_cheb2]=cheb2ord(Fp/(Fm/2),(Fp+BT)/(Fm/2),Rp,Rs);
[B_cheb2,A_cheb2]=cheby2(N_cheb2,Rp,Wn_cheb2);
tarea1(B_cheb2,A_cheb2,'Tchebyshev 2 filter');

[N_ellip,Wn_ellip]=ellipord(Fp/(Fm/2),(Fp+BT)/(Fm/2),Rp,Rs);
[B_ellip,A_ellip]=ellip(N_ellip,Rp,Rs,Wn_ellip);
tarea1(B_ellip,A_ellip,'Elliptic filter');

fvtool(B_butt,A_butt,B_cheb1,A_cheb1,B_cheb2,A_cheb2,B_ellip,A_ellip)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=tarea1(B,A,ttl)
figure

subplot(2,2,1)
impz(B,A);

[H,w] = freqz(B,A);
subplot(4,2,2)
plot(w/pi,abs(H))
subplot(4,2,4)
plot(w/pi,angle(H))

subplot(2,2,3)
zplane(B,A);

subplot(2,2,4)
grpdelay(B,A);
title(ttl,'Interpreter','latex');
end