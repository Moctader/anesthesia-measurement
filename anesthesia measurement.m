clc
clear all

load('521282S_data_2.mat')
load('521282S_data_2.mat', 'Fs')
load('521282S_data_2.mat', 'signal')
load('521282S_data_2.mat', 't')

%% 
Fs = 200;  % Sampling Frequency

N    = 800;      % Order
Fc1  = 1;        % First Cutoff Frequency
Fc2  = 4;        % Second Cutoff Frequency
flag = 'scale';  % Sampling Flag
% Create the window vector for the design algorithm.
win = hamming(N+1);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, [Fc1 Fc2]/(Fs/2), 'bandpass', win, flag);
Hd = dfilt.dffir(b);
delta_fil=filter(Hd,signal);
%% 
Fs = 200;  % Sampling Frequency

N    = 800;      % Order
Fc1  = 4;        % First Cutoff Frequency
Fc2  = 8;        % Second Cutoff Frequency
flag = 'scale';  % Sampling Flag
% Create the window vector for the design algorithm.
win = hamming(N+1);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, [Fc1 Fc2]/(Fs/2), 'bandpass', win, flag);
Ha = dfilt.dffir(b);
theta_fil=filter(Ha,delta_fil);

%% 

Fs = 200;  % Sampling Frequency

N    = 800;      % Order
Fc1  = 8;        % First Cutoff Frequency
Fc2  = 12;       % Second Cutoff Frequency
flag = 'scale';  % Sampling Flag
% Create the window vector for the design algorithm.
win = hamming(N+1);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, [Fc1 Fc2]/(Fs/2), 'bandpass', win, flag);
He = dfilt.dffir(b);
alpha_fil=filter(He,theta_fil);

%% 

Fs = 200;  % Sampling Frequency

N    = 800;      % Order
Fc1  = 12;       % First Cutoff Frequency
Fc2  = 25;       % Second Cutoff Frequency
flag = 'scale';  % Sampling Flag
% Create the window vector for the design algorithm.
win = hamming(N+1);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, [Fc1 Fc2]/(Fs/2), 'bandpass', win, flag);
Hf = dfilt.dffir(b);
beta_fil=filter(Hf,alpha_fil);
%%
figure()
subplot(5,1,1)
plot(t,signal);
subplot(5,1,2)
plot(t,delta_fil)
subplot(5,1,3)
plot(t,theta_fil)
subplot(5,1,4)
plot(t,alpha_fil)
subplot(5,1,5)
plot(t,beta_fil)

%%
[Env1, Env2]=envelope(signal,30*Fs,'rms');
[Env3, Env4]=envelope(delta_fil,30*Fs,'rms');
[Env44, Env5]=envelope(theta_fil,30*Fs,'rms');
[Env6, Env7]=envelope(alpha_fil,30*Fs,'rms');
[Env8, Env9]=envelope(beta_fil,30*Fs,'rms');

%%
figure()
subplot(5,1,1)
hold on
plot(t,signal);
plot(t,Env1,'r')
plot(t,Env2,'r')

subplot(5,1,2)
hold on
plot(t,delta_fil)
plot(t,Env3,'r')
plot(t,Env4,'r')

subplot(5,1,3)
hold on
plot(t,theta_fil)
plot(t,Env44,'r')
plot(t,Env5,'r')

subplot(5,1,4)
hold on
plot(t,alpha_fil)
plot(t,Env6,'r')
plot(t,Env7,'r')

subplot(5,1,5)
hold on
plot(t,beta_fil)
plot(t,Env8,'r')
plot(t,Env9,'r')


%% spectrogram and relative power
figure()
subplot(3,1,1)
[S,F,T,P] = spectrogram(signal,30*Fs,29*Fs,[0.1:0.1:32],Fs);

imagesc(T/60,F,log10(P),[-7 -3]);axis xy;

a=sum(P(10:40,:))./sum(P);
b=sum(P(40:80,:))./sum(P);
c=sum(P(80:120,:))./sum(P);
d=sum(P(120:250,:))./sum(P);

subplot(3,1,2)
hold on
plot(T/60,a,'r')
hold on
plot(T/60,b,'b')
hold on 
plot(T/60,c,'g')
hold on 
plot(T/60,d,'c')

%% spectral entrphy
SE=zeros(1,length(T))
for i=1:length(T)
    Pnorm=P(:,i)/sum(P(:,i));
    SE(i)=-sum(Pnorm.*log2(Pnorm))/log2(length(Pnorm));
end
subplot(3,1,3);
plot(T/60,SE)

















