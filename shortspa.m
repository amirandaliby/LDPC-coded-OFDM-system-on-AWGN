%LDPC coded OFDM system in AWGN
clc
clear all
M=2;%assigning the M-ary for psk modulation
cp=0; % Cylic prefix length

H=[1 1 0 1 0 0;0 1 1 0 1 0;0 0 1 1 0 1]; %parity chk matrix
G=[0 1 1 1 0 0;1 1 0 0 1 0;1 1 1 0 0 1]; %generator polynomials
%size
[a,b] = size(H);
[c,d] = size(G);
datasize=c;
Nfft=d;
%code rate, row weight and column weight
%Wr=3, %Wc=1-2, R=irregular;  
frame=100;

SNR1=0:1:20;%assigning the value of snr in decibels
SNR=10.^(SNR1/10);%converting the snr db value to dynamic range

%==== ITU Pedestrian channel  B ====
h=zeros(1,43);
p=[-3.92 -4.82 -8.82 -11.92 -11.72 -27.82];
cc=sqrt(10.^(p./10));
ind=[1 3 10 15 27 43];
h(ind)=cc;
KK=length(h);        

for i=1:length(SNR)
        error=0;
        for jj=1:frame;%MONTE CARLO LOOP

%---------------------------------TRANSMITTER----------------------------------  
            msg=randint(1,c,M); % Random data

%----------LDPC Encoding------------
            code=mod(msg*G,2);                    % modulo-2 Encoding
%----------interleaving----------
            intrlvcode=randintrlv(code,4831);
%-----------PSK modulation------------------
            modulator=modem.pskmod('M',2,'PhaseOffset',0,'SymbolOrder','gray','InputType','integer');%definng the 8-psk modulator objec
            mod_code=modulate(modulator,intrlvcode);%modulate with the input
%-----------IFFT & CP------------------
            y_ifft=ifft(mod_code,Nfft)*sqrt(Nfft);%performing the inverse fast fourier transform
            y_ifft_cp=[y_ifft(end-cp+1:end) y_ifft];%addition of cyclic prefix with last k sample

%---------------------------------CHANNEL----------------------------------  

 %----------mulipath fading----------  
            U=conv(h,y_ifft_cp);
%-------------AWGN Noise--------------------
            sigma=sqrt(1/(SNR(i)*log2(M)*2*0.5));%to compute the value of sigma 
            w=sigma*(randn(1,length(U)));%generating complex additive gaussian noise
           
            ytr=w+U;%adding the complex additive gaussian noise to the y_ifft_cp


%---------------------------------RECEIVER---------------------------------
%---------- Discard last Nch-1 received values resulting from convolution
            ytr(end-KK+2:end) = [];
%-----------IFFT & CP removal------------------
            ytr(1:cp) = []; %to remove the gaurd interval
            Y=fft(ytr,Nfft)/sqrt(Nfft);%performing the fast fourier transform
%-----------Zero Force (ZF) one-tap equalizer--------------------
            H=fft(h,Nfft);
            Yeq=Y./H;
%-----------PSK de-modulation------------------
%             demodulator=modem.pskdemod('M',2,'PhaseOffset',0,'SymbolOrder','gray');%defining the 8 psk demodulator object
%             demod=demodulate(demodulator,Y);%perform demodulation
%-----------deinterleaving--------------------
            deintrlvdemod=randdeintrlv(Yeq,4831); 
            
            LLR=(2*deintrlvdemod)/(sigma*sigma);
           
%-----------sum product algorithm Decoding----
%----Stage 1 : Initialization(Q=H*q)

            for x=1:a
             for y=1:b
                     Q(x,y)= H(x,y)*LLR(y);
             end
            end
            
            for iterations=1:40
 
%----Stage 2 : Horizontal step
                Tanh = tanh (Q./2); %Tanh of Q
                varb= Tanh; %defining a variable
                varb(varb==0)=1;
                cvec = prod (varb,2); % column vector (product of all elements in the same row)
                prct=repmat (cvec,1,b);  %repeating the column vector N
      

                for x=1:a
                    for y=1:b
                        if varb(x,y)~=1
                            L(x,y)= prct(x,y)/varb(x,y); %excluding same elemnt of the product
                           
                            if L(x,y)<-0.9999
                                L(x,y)=-0.9999;
                            elseif L(x,y)>0.9999
                                L(x,y)=0.9999;
                            end
                            
                            R(x,y)= log ((1+L(x,y))/(1-L(x,y)));   
                       
                        else
                            R(x,y)=0;
                        end
                    end 
                end
 
 
softout = zeros(1,b);
%-----------LDPC Hard decision-----------------
 
                vertical_sum=sum (R);
                for y=1:b
                    softout(y)= LLR(y)+ vertical_sum(y);
                    if softout(y) < 0  
                     out(y)= 1;
                    else
                     out(y)= 0;
                    end
                end
 
                
%-----------Break point-----------------               
                stop=mod(out*H',2);
                if stop == zeros(1,a)
                    break;
                end
 
%----stage 3 : Vertical step
                    for x=1:a
                     for y=1:b
                         if H(x,y)~=0
                            Q(x,y)=LLR(y)+ vertical_sum(y)-R(x,y); %updated Q
                         else
                            Q(x,y)=0;
                         end
                     end
                    end
                
 
    
            end
%--------------------------------------------------------------------------
[err ratio] = biterr(msg,out(4:6));%comparing the error in the data received to the input
error=error+err;
            
       end
error1(i)=error/(3*log2(M)*frame);
end


%computing BER performance in AWGN
p=0.5*erfc(sqrt(SNR));

%plotting the results
figure(1)
semilogy(SNR1,p);
hold on
semilogy(SNR1,error1,'ro-');
xlabel('SNR in dB')
ylabel('Bir error rate')
title('Bit error rate vs SNR for AWGN channel cp=64')
grid on
