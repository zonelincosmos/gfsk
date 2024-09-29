clc;clear all;close all;
%% GFSK Baseband Signal Settings 
% 20240929
symbolrate           = 1e6;                  
freq_deviation     = 160e3;   
modulation_idx   = (2*freq_deviation)/symbolrate; 

%% Gaussian Filter 
BT =  0.5;                   %  3-dB bandwidth-symbol time product
Span = 3;                       % Number of symbols
OversampleFactor = 8;     %  Samples per symbol period (oversampling factor)
Fs = symbolrate*OversampleFactor;

h = gaussdesign(BT, Span, OversampleFactor); 
% h = [1 zeros(1,48)]
% fvtool(h)

%% Transimitted data
SymbolCnt = 1202;
data = randi([0, 1],SymbolCnt,1)
% load data.m;   
% data = ones(120,1);
% data = zeros(120,1);
% data = [zeros(1,51) ones(1,1151)].';
% data(1:800) = 1;

%% Oversample 
b2= 2*data - 1;
b2(SymbolCnt,1) = b2(1,1);
bb2 = filter(ones(OversampleFactor,1), 1, upsample(b2,OversampleFactor));

%% Filtering
M = SymbolCnt*OversampleFactor; 
rr_g = ifft(fft(h.',M).*fft(bb2,M));

%% GFSK modulation
rr_m = zeros(1, M);
if(1)
    
    rr_m(1) = rr_g(1);
    for i= 2:1:M
        rr_m(i)=rr_m(i-1) + rr_g(i)/OversampleFactor;
    end
    fai=2*pi*(modulation_idx/2)*rr_m; 
    
else
    
    rr_m(1) = rr_g(1)*(1/Fs);
    for i= 2:1:M
        rr_m(i)=rr_m(i-1) + rr_g(i)*(1/Fs);
    end
    fai=2*pi*freq_deviation*rr_m; % 2*pi*f_deviation*t
    
end

DataComplex = cos(fai) + 1i*sin(fai);
%% Plot Waveform
subplot(3,1,1)
plot(real(DataComplex(1:1000)), 'LineWidth', 2)
hold on
plot(bb2(1:1000))
axis([0 1000 -1.5 1.5])
grid on

subplot(3,1,2)
plot(imag(DataComplex(1:1000)), 'LineWidth', 2)
hold on
plot(bb2(1:1000))
axis([0 1000 -1.5 1.5])
grid on
%% Plot Spectrum
FreqRange = [-Fs/2: Fs/M : Fs/2-Fs/M];
fftout = fftshift(fft(DataComplex));
subplot(3,1,3)
plot(FreqRange, 20*log10(abs(fftout)))
grid on
axis([-Fs/2 Fs/2 -20 80])


% sa1 = dsp.SpectrumAnalyzer('SampleRate', Fs);
% sa1.YLimits = [-10,90]
% sa1(DataComplex.')

%%
figure
plot(fai(1:1000))
hold on
plot(bb2(1:1000))
grid on
%% GFSK Demodulation
% Freq1 = +freq_deviation;
% Freq2 = -freq_deviation;
% time = [0: 1/Fs: (M - 1)*1/Fs];
% Freq1CosWaveform = cos(2*pi*Freq1*time);
% Freq1SinWaveform = sin(2*pi*Freq1*time);
% Freq2CosWaveform = cos(2*pi*Freq2*time);
% Freq2SinWaveform = sin(2*pi*Freq2*time);
% 
% q1 = Freq1CosWaveform.*real(DataComplex);
% q2 = Freq2CosWaveform.*real(DataComplex);
% plot(q1(1:1000))
% plot(abs(q1(1:1000)-1i*q2(1:1000)))
% hold on
% plot(bb2(1:1000))


