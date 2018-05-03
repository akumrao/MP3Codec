% Project - 2
% Pyschoacoustics 
% MP3 Codec

clear all
close all

[s,fs] = audioread('/root/Kashmir.wav');
output = zeros(length(s),1);
bitRate = 470000;
%output = zeros(1,length(s_Padded));

%% Absolute Threshold of Listening (dB-SPL)

f=1:fs/2;
tqf = 3.64*((f./1000).^-0.8) - 6.5*exp((-0.6)*(((f./1000)-3.3).^2))+10^-3*((f./1000).^4);

figure ;
plot(f, tqf);

title ("Absolute Threshold of Listening (dB-SPL)");
 
xlim ([20, 20000]);
ylim ([0, 100]);


semilogx(f,tqf)
 grid on;
axis([50 17000 -10 100])
title('Abolute threshold of hearing in quiet');
 xlabel('Frequency(Hz)');
ylabel('Sound Pressure Level,SPL(dB)');

%% Critical Bandwidth

% Critical Bandwidth
BWc = 25+75*(1+1.4*((f./1000).^2).^0.69);

figure ;
plot(f, BWc);

title ("Critical Bandwidth hz-hz)");
 
xlim ([20, 20000]);
ylim ([0, 6000]);

% Bark Scale Conversion
zfc = 13*atan(0.00076.*f)+3.5*atan((f/7500).^2);
% Centre Frequencies
f_centre = [50 150 250 350 450 570 700 840 1000 4800 5800 7000 8500 10500 13500 19500];
Bandwidths = [20 100 200 300 400 510 630 770 920 1080 1270 1480 1720 2000 2320 2700 3150 3700 4400 5300 6400 7700 9500 12000 15500 22050];
for i=1:1:25
    BarkBand(i) = Bandwidths(i+1)-Bandwidths(i);
end
 figure()
% semilogx(f,BWc);
 grid on;
 axis([0 22000 0 6000]);
 title('Critical BandWidth');
 xlabel('Frequency,(Hz)');
 ylabel('Critical Bandwatch(Hz)');
 figure()
 plot(f,zfc);
 title('Bark Scale Conversion');
xlabel('Frequency,(Hz)');
 ylabel('Critical Bandwidth(Bark)');
 grid on;
 axis([0 22000 0 25]);

%% Step -1 MPEG Codec
%Grabbing Frames
%Window size
N=512;
%Hanning Window
win=hanning(N);
%Shift by 128 samples 75% overlap
shift=N/16;
l = floor((length(s) - N)/shift);windSig = zeros(512,1);
nFrames=length(s)/(256)-1;

%%
for frames=1:nFrames
    startframe= (frames-1)*256+1;
    endframe=startframe+N-1;
    windSig = s(startframe:endframe).*win;
    
%% 100th Fame
%% figure()
%% plot(windSigcomplete(100,:));
%% axis tight;
%% grid on;
%% title('100th Frame');
%% Computing Pk
f=1:fs/(N):fs/2;
% PSD Estimate ( Power Spectral Density) 
pk = 90.302 + 20*log10((abs(fft(windSig))));
pk = pk(1:256,1);
tqf = 3.64*((f./1000).^-0.8) - 6.5*exp((-0.6)*(((f./1000)-3.3).^2))+10^-3*((f./1000).^4);
 figure()
 plot(f,pk,f,tqf);
 axis([90 20000 -5 130]);
 title('PSD estimate');
 xlabel('Frequency (Hz)');
 ylabel('SPL (dB)');
 hold on
 f1=1:fs/2;
 for b = 1:length(Bandwidths)
     stemx(Bandwidths(b)) = 200;
 end
 stem(f1,stemx,'k', 'LineStyle', ':', 'marker','none');
 hold on
%% set(gca,'xscal','log')
%% Bark Bands

%Centre Frequencies
f_centre = [50 150 250 350 450 570 700 840 1000 4800 5800 7000 8500 10500 13500 19500];
Bandwidths = [20 100 200 300 400 510 630 770 920 1080 1270 1480 1720 2000 2320 2700 3150 3700 4400 5300 6400 7700 9500 12000 15500 22050];
for i=1:1:25
    BarkBand(i) = Bandwidths(i+1)-Bandwidths(i);
end
% Bark Scale Conversion
zfc = 13*atan(0.00076.*f)+3.5*atan((f/7500).^2);
% Critical Bandwidth
BWc = 25+75*(1+1.4*((f./1000).^2).^0.69);
zfc_Bandwidths = 13*atan(0.00076.*Bandwidths)+3.5*atan((Bandwidths/7500).^2);
 figure()
plot(zfc,pk,zfc,tqf);
 axis([0.9 25 1 130]);
 title('PSD estimate - Barks');
 xlabel('Bark (Hz)');
 ylabel('SPL (dB)');
 hold on;
 fbark=1:25;
 for b = 1:length(fbark)
     stemx_bark(fbark) = 200;
 end
 
 stem(fbark,stemx_bark,'k', 'LineStyle', ':', 'marker','none');
 hold on;

%% Step -2 Identification of Tonal Maskers 

lmo=1;
for k=2:256
    if (2 < k && k < 63)
	      dk =  2;
	   elseif (63 <= k && k < 127)
	      dk = [2 3];
	   elseif (127 <= k && k < 250)
	      dk = [2 6];
    end
    if (2 < k && k < 63)
         if pk(k)>pk(k+1) && pk(k)>pk(k-1) && pk(k)>pk(k+dk)+7 && pk(k)>pk(k-dk)+7
        st(lmo)=pk(k);
        kval(lmo) = k; 
        lmo=lmo+1;
         end
    end
    if (63 <= k && k < 127)
    if pk(k)>pk(k+1) && pk(k)>pk(k-1) && pk(k)>pk(k+dk(1,1))+7 && pk(k)>pk(k-dk(1,1))+7 ...
           && pk(k)>pk(k+dk(1,2))+7 && pk(k)>pk(k-dk(1,2))+7 
        st(lmo)=pk(k);
        kval(lmo) = k;
        lmo=lmo+1;
    end
    end
       if (127 <= k && k < 249)
       if pk(k)>pk(k+1) && pk(k)>pk(k-1) && pk(k)>pk(k+dk(1,1))+7 && ...
               pk(k)>pk(k-dk(1,1))+7 && pk(k)>pk(k+dk(1,2))+7 && ...
               pk(k)>pk(k-dk(1,2))+7 && pk(k)>pk(k+(dk(1,2)-1))+7 && ...
               pk(k)>pk(k-(dk(1,2)-1))+7 && pk(k)>pk(k+(dk(1,2)-2))+7 && ...
               pk(k)>pk(k-(dk(1,2)-2))+7 && pk(k)>pk(k+(dk(1,2)-3))+7 && ...
               pk(k)>pk(k-(dk(1,2)-3))+7 
            st(lmo)=pk(k);
            kval(lmo) = k;
            lmo=lmo+1;
        end
        end
end
%% Plotting Raw Tonal Maskers - Just peaks, adjacent ones aren't summed ( not the PTM!! This is ST! )

knew=[];
sd =1;
for k=1:256
    if sd>length(st)
        break;
    elseif k==kval(sd)
        knew(k)=st(sd);
        sd=sd+1;
    end
end
for k=length(knew)+1:256
    knew(k) = 0;
end
 figure()
 plot(zfc,pk,'r', zfc,tqf,'-',zfc,knew,'x');
 axis([0.9 25 -5 130]);
 title('Spectral Peaks');
 xlabel('Bark (Hz)');
 ylabel('SPL (dB)');
 hold on;
 
 stem(fbark,stemx_bark,'k', 'LineStyle', ':', 'marker','none');
 hold on
%% Computing PTM and placing it at the exact k value

for hk=1:length(kval)
ptm(kval(1,hk)) = 10*log10((10^(0.1*pk((kval(1,hk)-1)))+10^(0.1*pk(kval(1,hk)))+10^(0.1*pk(kval(1,hk)+1))));
end
%% Plotting PTM 

for k=length(ptm)+1:256
    ptm(k) = 0;
end

 figure()
 plot(zfc,pk,zfc,tqf,'-',zfc,ptm,'x');
 axis([0.9 25 1 130]);
 title('Tonal Maskers');
 xlabel('Bark (Hz)');
 ylabel('SPL (dB)');
 hold on;
 stem(fbark,stemx_bark,'k', 'LineStyle', ':', 'marker','none');
% hold on

%% Computing Geometric mean kbar
 
counter = 1;
count = 1;
incr = round((fs)/(N));
for hq=1:length(Bandwidths)-1
    l=Bandwidths(hq);
    u=Bandwidths(hq+1);
    ilp=1;
    kbar =1;
    while (incr*counter>=l && incr*counter<=u)
       value = incr*counter;
        kbar = incr*counter*kbar;
       counter = counter+1;
       ilp = ilp+1;
    end
    kfu(count) = kbar;
    kbar_array(count) = round(nthroot(kbar,ilp-1));
    count = count+1;
end

kbar_array(1,25) = 0;

%% Calculating PNM 

crap = 256/25;
count = 1;
for val= 1:round(crap):256
   new(count) = val; 
   count=count+1;
end
count_final=1;
for val=1:length(new)-1
   lnew=new(val);
    unew=new(val+1);
    pk_noisemaskers = 0;
    for counter = lnew:unew
        pk_noisemaskers = (10.^(0.1*(pk(counter,1)))+pk_noisemaskers);
    end
    pk_noisemaskers_final(count_final) = pk_noisemaskers;
    count_final = count_final+1;
end

%%

kbar_array_krounded = round(kbar_array/86);
count=1;
for fgh = 1:255
    if kbar_array_krounded(1,count) == fgh
        pnm(fgh) = 10*log10(pk_noisemaskers_final(1,count));
        count = count+1;
    end  
end

%%

for k=length(ptm)+1:256
    ptm(k) = 0;
end

for k=length(pnm)+1:256
    pnm(k) = 0;
end
% figure()
% plot(zfc,pk,zfc,tqf,'-',zfc,ptm,'x');
% axis([0.9 25 -1 130]);
% title('Both Tonal and Noise Maskers');
% xlabel('Bark (Hz)');
% ylabel('SPL (dB)');
% hold on;
% 
% stem(fbark,stemx_bark,'k', 'LineStyle', ':', 'marker','none');
% hold on
% 
% plot(zfc,pnm,'o');
% hold on

%% Checking for components from k which can be removed

lk=1;
for k=2:255
    if (2 < k & k < 63)
	      dk =  2;
	   elseif (63 <= k & k < 127)
	      dk = [2 3];
	   elseif (127 <= k & k < 250)
	      dk = [2 6];
    end

    if (2 < k & k < 63)
         if pk(k)>pk(k+1) && pk(k)>pk(k-1) && pk(k)>pk(k+dk)+7 && pk(k)>pk(k-dk)+7
            noiseMaskKVal(lk) = k-1;
            noiseMaskKVal(lk+1) = k;
            noiseMaskKVal(lk+2) = k+1;
            lk=lk+3;
         end
    end
    if (63 <= k & k < 127)
    if pk(k)>pk(k+1) && pk(k)>pk(k-1) && pk(k)>pk(k+dk(1,1))+7 && pk(k)>pk(k-dk(1,1))+7 ...
           && pk(k)>pk(k+dk(1,2))+7 && pk(k)>pk(k-dk(1,2))+7 
            noiseMaskKVal(lk) = k-3;
            noiseMaskKVal(lk+1) = k-2;
            noiseMaskKVal(lk+2) = k-1;
            noiseMaskKVal(lk+3) = k;
            noiseMaskKVal(lk+4) = k+1;
            noiseMaskKVal(lk+5) = k+2;
            noiseMaskKVal(lk+6) = k+3;
            lk=lk+7;    
    end
    end
       if (127 <= k & k < 249)
       if pk(k)>pk(k+1) && pk(k)>pk(k-1) && pk(k)>pk(k+dk(1,1))+7 && ...
               pk(k)>pk(k-dk(1,1))+7 && pk(k)>pk(k+dk(1,2))+7 && ...
               pk(k)>pk(k-dk(1,2))+7 && pk(k)>pk(k+(dk(1,2)-1))+7 && ...
               pk(k)>pk(k-(dk(1,2)-1))+7 && pk(k)>pk(k+(dk(1,2)-2))+7 && ...
               pk(k)>pk(k-(dk(1,2)-2))+7 && pk(k)>pk(k+(dk(1,2)-3))+7 && ...
               pk(k)>pk(k-(dk(1,2)-3))+7 
            noiseMaskKVal(lk) = k-6;
            noiseMaskKVal(lk+1) = k-5;
            noiseMaskKVal(lk+2) = k-4;
            noiseMaskKVal(lk+3) = k-3;
            noiseMaskKVal(lk+4) = k-2;
            noiseMaskKVal(lk+5) = k-1;
            noiseMaskKVal(lk+6) = k;
            noiseMaskKVal(lk+7) = k+1;
            noiseMaskKVal(lk+8) = k+2;
            noiseMaskKVal(lk+9) = k+3;
            noiseMaskKVal(lk+10) = k+4;
            noiseMaskKVal(lk+11) = k+5;
            noiseMaskKVal(lk+11) = k+6;
            lk=lk+13;
        end
        end
end

%% Elimination of noise maskers which are tonal as well

% if kbarrounded == noisemaskKvalues.. remove them 
count = 1;
[Lia,Locb] = ismember(kbar_array_krounded,noiseMaskKVal);

for m = 1:256
    if pnm(m) ~= 0 
        pnm_final(count) = pnm(m);
        count = count+1;
    end
end
pnm_final(25) = 0;

%%
new1= [];
count =1;
for m = 1:length(Lia)
    if Lia(m) == 1
        new1(count) = m;
        count = count +1;
    end
end
for j = 1:length(new1)   
    pnm_final(new1(j)) = 0;
end

%% Re sizing the pnm_final

count =1;
for m = 1:length(pnm)
   
    if pnm(m) ~= 0
        new12(count) = m;
        count = count +1;
    end
    
end
new(25) = 256;
new_array = zeros(1,256);
for j = 1:length(new12)   
    new_array(new12(j)) = pnm_final(j);
end

%%

for k=length(ptm)+1:256
    ptm(k) = 0;
end

% figure()
% plot(zfc,pk,zfc,tqf,'-',zfc,ptm,'x');
% 
% title('Elimination Noise maskers and also anything below absolute threshold');
% xlabel('Bark (Hz)');
% axis([0.9 25 -1 130]);
% ylabel('SPL (dB)');
% hold on;
% 
% stem(fbark,stemx_bark,'k', 'LineStyle', ':', 'marker','none');
% hold on
% %
% plot(zfc,new_array,'o');
% hold on

%% Step 3 - Decimation 

for m = 1:256
    
    if ptm(m)<tqf(m)
        ptm(m) = 0;
    end
end

for m = 1:256
    
    if new_array(m)<tqf(m)
        new_array(m) = 0;
    end
end

%% Array containing both PTM and PNMs

% ptm - Tonal Maskers
% new_array - Noise Maskers

appended_array = zeros ( 1,256);

for m = 1:256
   
    if ptm(m)>0
        appended_array(m) = ptm(m) ;
    end
    if new_array(m)>0
        appended_array(m) = new_array(m); 
    end
end

%% Sliding Bark Window - PTM

[peaks,locations]=findpeaks(appended_array,'MinPeakDistance',4);

% Making array of equal size
new_array_slidingwind = zeros(1,256);
for j = 1:length(locations)   
    new_array_slidingwind(locations(j)) = peaks(j);
end

%% Getting back tonal and noise maskers as they were 

for m=1:256
    if ptm(m) == new_array_slidingwind(m)
       ptm_slide_wind(m) = new_array_slidingwind(m); 
    else 
        ptm_slide_wind(m) = 0;
    end 
end

for m=1:256
   
    if new_array(m) == new_array_slidingwind(m)
       pnm_slide_wind(m) = new_array_slidingwind(m); 
    else 
        pnm_slide_wind(m) = 0;
    end
    
end

%Plotting
% 
% figure()
% plot(zfc,pk,zfc,tqf,'-',zfc,ptm_slide_wind,'x');
% 
% title('Sliding Bark Window to Eliminate tonal or noise maskers');
% xlabel('Bark (Hz)');
% axis([0.9 25 -1 130]);
% ylabel('SPL (dB)');
% hold on;
% 
% stem(fbark,stemx_bark,'k', 'LineStyle', ':', 'marker','none');
% hold on
% 
% plot(zfc,pnm_slide_wind,'o');
% hold on

%plot(ttm_new)

%% Convertion of ks to i's

 i=[];
 for k = 1:1:232
     if k>=1 && k<= 48
           i(k)=k;
    elseif k>=49 && k<=96
        i(k) = k+(rem(k,2));
    elseif k>=97 && k<=232
        i(k) = k + 3 - (rem((k-1),4));
    end
 end

for m = length(i)+1:256
   i(m) = m; 
end

%% Decimation of ptm and pnms to i's

ptm_iterms = zeros(1,256);
for j = 1:length(i)   
    ptm_iterms(j) = ptm_slide_wind(i(j));
end

pnm_iterms = zeros(1,256);
for j = 1:length(i)   
    pnm_iterms(j) = pnm_slide_wind(i(j));
end
%% Normalizing to z's 
kbar_barks = 13*atan(0.00076.*kbar_array)+3.5*atan((kbar_array/7500).^2);
i_freq = (fs/N)*i;
zi = 13*atan(0.00076.*i_freq)+3.5*atan((i_freq/7500).^2);
locations_bark = locations/10.24;

%% Spreading Functions SF(i,j)

for n = 1:length(locations)
for m= 1:length(zi)
     dz(n,m) = zi(m) - locations_bark(1,n);
    if -3<=dz(n,m) && dz(n,m)<-1
        sfij(m,n) = 17*dz(n,m)-(0.4*ptm_iterms(1,i(1,locations(1,n))))+11;
    elseif -1<=dz(n,m) && dz(n,m)<0
        sfij(m,n) = (0.4*(ptm_iterms(1,i(1,locations(1,n))))+6)*dz(n,m);
    elseif 0<=dz(n,m) && dz(n,m)<1
        sfij(m,n) = -17*dz(n,m);
    elseif 1<=dz(n,m) && dz(n,m)<8
        sfij(m,n) = (0.15*ptm_iterms(1,i(1,locations(1,n)))-17)*dz(n,m)-0.15*(ptm_iterms(1,i(1,locations(1,n))));
    end
end
end
sfij(sfij==0) = NaN;
[l,g] = size(sfij);

%% Calculating locations of the final tonal maskers and noise maskers
count =1;
for m = 1:256
   if ptm_slide_wind(1,m)>0
       tm_locations(count) = m;
       count = count+1;
   end  
end

count = 1;
for m = 1:256
   if pnm_slide_wind(1,m)>0
       nm_locations(count) = m;
       count = count+1;
   end  
end

%% Threshold Tonal Maskers

tm_locations_freq = 86*tm_locations;
tm_locations_bark=13*atan(0.00076.*tm_locations_freq)+3.5*atan((tm_locations_freq/7500).^2);

for n = 1:length(tm_locations)
for m= 1:l 
    ttm(m,n) = ptm_slide_wind(1,tm_locations(1,n))-(0.275*tm_locations_bark(1,n))+sfij(m,n)-6.025;
end
end

%% Resizing Tonal Thresholds to get the plot

[maxval_ttm,maxval_location_ttm] = max(ttm);
ttm_new = zeros(size(ttm'));
[m,n] = size(ttm_new);
for g = 1:n
for l = 1:m
   if g-(maxval_location_ttm(1,l)-tm_locations(1,l))-1<0
       ttm_new(g,l) = 0;
   else
       ttm_new(g-(maxval_location_ttm(1,l)-tm_locations(1,l)),l)=ttm(g,l);
   end
end
end
ttm_new(ttm_new==0) = NaN;
nm = length(tm_locations);
ttm_final = zeros(256,nm);
[m,n] = size(ttm_new);
for l = 1:n
    for g=1:m
        ttm_final(g,l) = ttm_new(g,l);
        
    end
end
ttm_final = ttm_final(1:256,1:nm);

%% Plotting Thresholds of Tonal Maskers

% figure()
% plot(zfc,pk,zfc,tqf,'-');
% 
% title('Tonal Masker Thresholds');
% xlabel('Bark (Hz)');
% axis([0.9 25 -1 110]);
% ylabel('SPL (dB)');
% hold on;
% 
% stem(fbark,stemx_bark,'k', 'LineStyle', ':', 'marker','none');
% hold on
% %
% plot(zfc,ptm_slide_wind,'x');
% hold on
% 
% plot(zfc,ttm_final);

%% Threshold Noise Maskers

[l,g] = size(sfij);
nm_locations_freq = 86*nm_locations;
nm_locations_bark=13*atan(0.00076.*nm_locations_freq)+3.5*atan((nm_locations_freq/7500).^2);

for n = 1:length(nm_locations)
for m= 1:l 
    tnm(m,n) = pnm_slide_wind(1,nm_locations(1,n))-(0.175*nm_locations_bark(1,n))+sfij(m,n)-2.025;
end
end

%% Resizing Noise Maskers

[maxval_tnm,maxval_location_tnm] = max(tnm);
tnm_new = zeros(size(tnm'));
[m,n] = size(tnm_new);
for g = 1:n
for l = 1:m
   if g-(maxval_location_tnm(1,l)-nm_locations(1,l))-1<0
       tnm_new(g,l) = 0;
   else
       tnm_new(g-(maxval_location_tnm(1,l)-nm_locations(1,l)),l)=tnm(g,l);
   end
end
end
tnm_new(tnm_new==0) = NaN;
nm = length(nm_locations);
tnm_final = zeros(256,nm);
[m,n] = size(tnm_new);
for l = 1:n
    for g=1:m
        tnm_final(g,l) = tnm_new(g,l);
        
    end
end
tnm_final = tnm_final(1:256,1:nm);
%% Plotiing Thresholds of Noise Maskers
% figure()
% plot(zfc,pk,zfc,tqf,'-');
% 
% title('Noise Masker Thresholds');
% xlabel('Bark (Hz)');
% axis([0.9 25 -1 110]);
% ylabel('SPL (dB)');
% hold on;
% 
% stem(fbark,stemx_bark,'k', 'LineStyle', ':', 'marker','none');
% hold on
% %
% plot(zfc,pnm_slide_wind,'o');
% hold on
% 
% plot(zfc,tnm_final);

%% Plotting thresholds on Barks
% figure()
% plot(zfc,pk,zfc,tqf,'-',zfc,ptm_slide_wind,'x');
% 
% title('Threholds of both Tonal and Noise Maskers');
% xlabel('Bark (Hz)');
% axis([0.9 25 -1 110]);
% ylabel('SPL (dB)');
% hold on;
% 
% stem(fbark,stemx_bark,'k', 'LineStyle', ':', 'marker','none');
% hold on
% 
% plot(zfc,pnm_slide_wind,'o',zfc,ptm_slide_wind,'x');
% hold on
% 
% plot(zfc,ttm_final,zfc,tnm_final);
% hold on

%% Computing Global Thresholds 

ttm_new(isnan(ttm_new))=0;
tnm_new(isnan(tnm_new))=0;
ttm_sum = sum(10.^(0.1*(ttm_new')));
tnm_sum = sum(10.^(0.1*(tnm_new')));
ttm_sum_final = zeros(1,256);
for l = 1:length(ttm_sum)
   ttm_sum_final(1,l) = ttm_sum(1,l); 
end
tnm_sum_final = zeros(1,256);
for l = 1:length(tnm_sum)
   tnm_sum_final(1,l) = tnm_sum(1,l); 
end

%%

for df = 1:length(tqf)
    tg(df) = 10*log10(10.^(0.1*tqf(df)) + tnm_sum_final(df) + ttm_sum_final(df));
end

%%
% figure()
% plot(zfc,pk,zfc,tqf,'-',zfc,ptm_slide_wind,'x');
% 
% title('GLOBAL THRESHOLD');
% xlabel('Bark (Hz)');
% axis([0.9 25 -1 110]);
% ylabel('SPL (dB)');
% hold on;
% 
% stem(fbark,stemx_bark,'k', 'LineStyle', ':', 'marker','none');
% hold on
% 
% plot(zfc,pnm_slide_wind,'o',zfc,ptm_slide_wind,'x');
% hold on
% 
% plot(zfc,ttm_final,zfc,tnm_final);
% hold on
% 
% plot(zfc,tg,'k');
% hold on

%% Bit Allocation

% bit allocation
    desSMNR=tg;
    actualBits=zeros(256,1);
    count=0;
    rounding_bitrate=round(bitRate/2*N/fs);
    while(count~=rounding_bitrate)
        [yn, ind]=min(desSMNR);
        actualBits(ind)=actualBits(ind)+1;
        desSMNR(ind)=desSMNR(ind)+6.02; % 1bit = 6.02dB
        count=count+length(ind);
    end
    dct_sn=dct(s(startframe:endframe).*win);
    bitCrushedDct=zeros(256,1);
    dct_snNN=zeros(256,1); %non-negative
    for lmo=1:256        
        if(dct_sn(lmo)>=0)
            dct_snNN(lmo)=1;
        end
    end
    dct_sn=abs(dct_sn);
    for lmo=1:256
        bit2Num=2^(actualBits(lmo)-1);
        bitCrushedDct(lmo)=round(bit2Num*dct_sn(lmo)/16)*16/bit2Num;
        if(dct_snNN(lmo)==0)
            bitCrushedDct(lmo)=-bitCrushedDct(lmo);
        end
    end
    idct_sn=idct(bitCrushedDct,N);
    tempFrame=zeros(length(s),1);
    tempFrame(startframe:endframe)=idct_sn;
    output=output + tempFrame;
end

soundsc(output,fs);

wavwrite(output,fs,'Reconstructed-Kashmir320')