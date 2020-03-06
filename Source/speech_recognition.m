clc;
clear;
close all;

InputDir = "Input";
Subfolders = dir(InputDir);
len = 1000000000;
thres= 0.5;
thres1 = 0.35;
OutputDir = "Output_recognition";

[test_rec,Fs]=audioread("Test/Test1.wav");
test_rec=test_rec';
sztest=size(test_rec);
if sztest(1)>1
    test_rec = (test_rec(1,:)+test_rec(2,:))/2;
end
% fft_test=genDFT(test_rec,Fs,length(test_rec))
% fft_test=fft_test/max(fft_test)
% plot(fft_test)
test_rec=test_rec/max(test_rec);
% len=min(len,length(test_rec))
% [test_power,test_f]=periodogram(test_rec,[],[],Fs,'power');
% test_power=test_power/max(test_power)
% maxm_corr=[]
% maxm_corr_freq=[]
avg_cohr=[];
for k = 3:length(Subfolders)
    Person = Subfolders(k).name;
    Files = dir(fullfile(InputDir + "/" + Person, "*.wav"));
    
    figure1 = figure('visible','off', 'Position', get(0, 'Screensize'));
    title(sprintf("Similarity with voice of %s", Person));    
    for j = 1:length(Files)
        [x, Fs] = audioread(InputDir + "/" + Person + "/" + Files(j).name);
        x = x';
        szx = size(x);
        if szx(1)>1
            x = (x(1, :)+x(2, :))/2;
        end
        x=alignsignals(x,test_rec);
        len = min(len, length(x));
        x=x/max(x);
%         fft_x=genDFT(x,Fs,length(x))
%         fft_x=fft_x/max(fft_x)
%         [x_power,x_f]=periodogram(x,[],[],Fs,'power');
%         x_power=x_power/max(x_power)
%         [corr_sim,lag]=xcorr(test_rec,x);
%         maxm_corr=[maxm_corr;max(corr_sim)]    
%  +       
%         [corr_sim_freq,lag_freq]=xcorr(test_power,x_power);
%         maxm_corr_freq=[maxm_corr_freq;max(corr_sim_freq)]         
%         [corr_fft,lag]=xcorr(fft_x,fft_test)
        
        Timing = Files(j).name;
        Timing = Timing(1:end-4);        
        ax(j)=subplot(2, 3, j);
%         hold on
%         plot(x_f,x_power)
%         plot(test_f,test_power,'r')
%         hold off
        if length(test_rec)<length(x)
            test_rec_1=[test_rec zeros(1,length(x)-length(test_rec))];
        elseif length(test_rec)>length(x)
            x=[x zeros(1,length(test_rec)-length(x))];
            test_rec_1=test_rec ;
        else 
            test_rec_1=test_rec;
        end
        [cross_spect_density,f]  = cpsd(x,test_rec_1,[],[],[],Fs);
        cross_spect_density=abs(cross_spect_density);
        cross_spect_density=cross_spect_density/max(cross_spect_density);
        [sx,f]=cpsd(x,x,[],[],[],Fs);
        sx=abs(sx);
        sx=sx/max(sx);
        [stest,f]=cpsd(test_rec_1,test_rec_1,[],[],[],Fs);
        stest=abs(stest);
        stest=stest/max(stest);
        
        f1=[];
        for i=1:length(f)
            if stest(i)>=thres
                f1=[f1;f(i)];
            end
        end
        
        [spect_coh,f2] = mscohere(x,test_rec_1,[],[],f1,Fs);
        spect_coh_new=[];
        ctr=0;
        for i=1:length(f)
            if stest(i)>=thres
                ctr=ctr+1;
                spect_coh_new=[spect_coh_new;spect_coh(ctr)];
            else
                spect_coh_new=[spect_coh_new;0];
            end
        end        
        avg_cohr=[avg_cohr;mean(spect_coh)];
        
%         hold on
%         plot(x)
%         plot(test_rec_1,'r')
%         hold off
%         ylabel('Amplitude')
%         xlabel('Time')
        hold on;
        plot(f,cross_spect_density);
        stem(f,spect_coh_new,'r');
        hold off;
        ylabel('Ampl (red is spectral coherence)');
        xlabel('Freq');
        xlim([0,1000])
        grid on
        title(sprintf("In the %s", Timing));

        fprintf(Person + " " + Files(j).name + " Processed\n");
    end
    FileDir = OutputDir + "/" + Person;
    if ~exist(FileDir, 'dir')
       mkdir(FileDir)
    end
    saveas(figure1, FileDir + "/" + "Cross Spectral Density and auto coherence", "png");    
end
avg_cohr(1:6)=sort(avg_cohr(1:6),'descend');
avg_cohr(7:12)=sort(avg_cohr(7:12),'descend');
avg_cohr(13:18)=sort(avg_cohr(13:18),'descend');

mean1=mean(avg_cohr(1:4));
mean2=mean(avg_cohr(7:10));
mean3=mean(avg_cohr(13:16));
max_mean=max(mean1,max(mean2,mean3));
disp(" ")
disp("The Voice matches with")
if  max_mean<=thres1
    disp("No One")
elseif mean1>=mean2 & mean1>=mean3
    disp("Mayank")
elseif mean2>=mean1 & mean2>=mean3
    disp("Poojitha")
elseif mean3>=mean2 & mean3>=mean2
    disp("Prerna")
end        
    

function y=genDFT(x, Fs, N)

Fx = fft(x, N);
m = abs(Fx);
shifted_m = fftshift(m);
y=20*log10(shifted_m/N);

end