prompt = 'Imput the Samples Per Bit ';
SamplesPerBit = input(prompt);

prompt2 = 'Imput the number of bits to Tx (will be taken from an array already defined)';
BitsForTx = input(prompt2);
%SamplesPerBit=100; %Input Samples per bit
%BitsForTx= 10;  %Number of bits fot Tx, this number represents the number of bits taked from an array of random data already defined in d;
%If desired, comment this line and add a raw input of bits in vector 'd'.
%It's started with 1,0... with the purpose of being illustrative

d=zeros(1,BitsForTx);
d1=[1 0 1 0 1 0 0 1 1 1 1 1 1 0 0 1 0 1 0 0 0 0 0 0]; % Bits sequence, if desired define this as a random vector. 
for v=1:1:BitsForTx
   d(v)=d1(v); 
end    
nbits = length(d);
nbitsarray = (0:1:nbits-1);
b=2*d-1; % Convert unipolar to bipolar
T=1; % Bit duration
Eb=T/2; % This will result in unit amplitude waveforms
fc=4/T; % Carrier frequency, numerator = number of periods per bir
t=linspace(0,nbits,SamplesPerBit*nbits); %Array of samples in time
N=length(t); % Number of total samples
Nsb=N/length(d); % Number of samples per bit, its computed again to ensure an INT number of waves
ff0=repmat(d',1,Nsb); % replicate each bit NumberOfSamples times
fb=repmat(b',1,Nsb); 
dataf1=ff0'; 
dataf1=dataf1(:)'; % Transpose matrix
bw=fb';
bw=bw(:)'; % Data sequence samples


w=sqrt(2*Eb/T)*cos(2*pi*fc*t); % carrier signal
bpsk_w=bw.*w; % modulated waveform


%PLOT
sgtitle('Modulación BPSK');
subplot(4,1,1);
stem(nbitsarray,d,'filled', 'linewidth',2); axis([0 nbits 0 1]) 
xlabel('Bits');
ylabel('Amplitd ');
grid on;


subplot(4,1,2);
plot(t,bw); axis([0 nbits -1.5 1.5]) 
xlabel('Tiempo [s]');
ylabel('Amplitud');
grid on;

subplot(4,1,3);
plot(t,w); %axis([0 5 -1.5 1.5])
xlabel('Tiempo [s]');
ylabel('Amplitud');
grid on;

subplot(4,1,4);
plot(t,bpsk_w,'-'); %axis([0 5 -1.5 1.5])
xlabel('Tiempo [s]');
ylabel('Amplitud');
grid on;

%%%DEMODULACION
y=0;
SARRAY=zeros(1,SamplesPerBit);

bwRX =  bpsk_w.*w; %Product of Received signal by Carrier W

b=SamplesPerBit;
c=1;
k=[];

for m=1:nbits %
y=bwRX(c:b); % Divide signal in number of samples per bit
q=cos(2*pi*((c:b)-1)*fc/Nsb); %Signal for Quadrature for sync Demod
i=-cos(2*pi*((c:b)-1)*fc/Nsb); %InPhase signal for sync Demod
a=y.*q; % received signal by Quarature, 
d=y.*i; % received signal by inPhase 
%These signals are updated down
addAppend= sum(a)-sum(d); %Two symbols
if addAppend>0
    p=ones(1,SamplesPerBit); %Decisión to all samples to 1
else
    p=zeros(1,SamplesPerBit); %Decision to all samples to 0
end

k=[k p]; %append values in k
b=b+SamplesPerBit; %nex value
c=c+SamplesPerBit; %next value
end

%PLOT

figure(2);
sgtitle('Demodulación BPSK');
subplot(4,1,1);
plot(t,w); 
xlabel('Tiempo [s]');
ylabel('Amplitud');
grid on;

subplot(4,1,2);
plot(t,bpsk_w,'-'); 
xlabel('Tiempo [s]');
ylabel('Amplitud');
grid on;

subplot(4,1,3); plot(bwRX);grid on; 
axis([0 nbits*SamplesPerBit -1.5 1.5]); 
xlabel('Muestras');
ylabel('Amplitud');
grid on;
subplot(4,1,4);plot(k,'LineWidth',1.5);grid on;
axis([0 nbits*SamplesPerBit -0.5 1.5]);
xlabel('Muestras');
ylabel('Amplitud');
grid on;
