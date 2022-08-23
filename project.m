clc;
clear all; %#ok<CLALL>
close all;

VAL1 = 1;
transmitted_msg='PAPRI-MAISHA-TANJIM-PRANTO-NUBAH';
disp(transmitted_msg)


% 1))

ascii=double(transmitted_msg); % Text to ASCII (decimal)
disp(ascii);


% 2))

bin=dec2bin(ascii); % Converting Information Message to bit
nibble=dec2nibble(ascii); % Grouped to 4-bit
bp=.000001;
% bit period
disp('Information at Transmitter : ');
disp(nibble');

% 3))

bit=[];
for n=1:1:length(bin)
    if bin(n)==1
        se=5*ones(1,100);
    elseif bin(n)==0
        se=zeros(1,100);
    end
    bit=[bit se]; %#ok<AGROW>
end
t1=bp/100:bp/100:100*length(bin)*(bp/100);
figure(1);
subplot(5,1,1);
plot(t1,bit,'lineWidth',2.5);
grid on;
axis([ 0 bp*length(bin) -.5 6]);
ylabel('amplitude(volt)');
xlabel(' time(sec)');
title('Message Signal (Binary)');


% 4))
% 16-QAM modulation
%Message Signal
N = length(nibble);
M = 16; %Modulation Order 16
data = nibble; %Message Infromation
s = qammod(data,M); %Modulated data

bit_rate = 10.^3;
f = bit_rate; %minimum carrier frequency
Tb = 1/bit_rate ; %bit duration
t = 0:(Tb/1000):Tb ;

%Transmitted Signal waveform
TxSig = [];
for l=1:length(data)
    Tx = real(s(l))*cos(2*pi*f*t) - imag(s(l))*sin(2*pi*f*t);
    TxSig = [TxSig Tx]; %#ok<AGROW>
end

%Wave forms of the Signal
subplot(5,1,2);
stairs(data,'lineWidth',2.5);
grid minor;
ylim([-0.5,M-0.5]);
xlim([0,N + 1]);
title('Message Signal (Nibbles - 4-bit)');
subplot(5,1,3);
plot(TxSig,'lineWidth',2.5);
grid minor;
title('16-QAM Modulated Signal');
xlim([0,N*10^3+N]);


% 5))

disp('*** AWGN Added ***');

%Channel Noise%
SNR = 30 + VAL1; %SNR in dB
r = awgn(s,SNR,'measured');

%Received Signal waveform
RxSig = awgn(TxSig,SNR,'measured'); %AWGN for Tx Waveforms

subplot(5,1,4);
plot(RxSig,'lineWidth',2.5);
grid minor;
title('16-QAM Modulated Signal with AWGN');
xlim([0,N*10^3+N]);

% Constellation Diagram of the Rx
scatterplot(r);
grid minor;
title('Constellation Diagram of 16-QAM')


% 6))

demodulated = qamdemod(r,16);
disp('Information at Reciver :');
disp(demodulated');


% 7))

recieved_signal = nibble2dec(demodulated)';
recieved_msg = char(recieved_signal);
disp("recieved message : " + recieved_msg);

figure(1);
subplot(5,1,5);
stairs(demodulated,'lineWidth',2.5);
grid minor;
ylim([-0.5,M-0.5]);
xlim([0,N + 1]);
title('Message Signal (Nibbles - 4-bit)');


% 8))

if strcmp(recieved_msg, transmitted_msg)
    disp("Recieved Message is Identical to Transmitted Message!")
else
    disp("Recieved Message didn't match with Transmitted Message!")
end

% 9))

%Channel Noise%
SNR = 12; %SNR in dB
r = awgn(s,SNR,'measured');

%Received Signal waveform
RxSig = awgn(TxSig,SNR,'measured'); %AWGN for Tx Waveforms

% Constellation Diagram of the Rx
scatterplot(r);
grid minor;
title('Constellation Diagram of 16-QAM')

demodulated = qamdemod(r,16);

recieved_signal = nibble2dec(demodulated)';
recieved_msg = char(recieved_signal);
disp("recieved message : " + recieved_msg);

if strcmp(recieved_msg, transmitted_msg)
    disp("Recieved Message is Identical to Transmitted Message!")
else
    disp("Recieved Message didn't match with Transmitted Message!")
end


%%% Functions %%%


function dn = dec2bin(dec)
p2=2.^(0:-1:-7); % 2^0,2^-1,...,2^-7
B=mod(floor(p2'*dec),2);
dn=reshape(B,1,numel(B));
end

function nibbles = dec2nibble(decimals)
L = length(decimals);
nibbles2 = zeros(L,2);
LSN = 15;   %LSN, ie 0000 1111
for n = 1 : L
    nibbles2(n,1) = bitshift(decimals(n),-4);
    nibbles2(n,2) = bitand(decimals(n),LSN);
end
nibbles = reshape(nibbles2.', [], 1);
end

function decimals = nibble2dec(nibbles)
tran = (reshape(nibbles, 2, [])).';
L = length(tran);
decimals = zeros(L,1);

for i = 1 : L
    MSB = tran(i,1)*2^4;
    decimals(i) = MSB+tran(i,2);
end
end
