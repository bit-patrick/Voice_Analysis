clc;
clear;
close all;

InputDir = "Input";
OutputDir = "Output";
Subfolders = dir(InputDir);
len = 1000000000;

for k = 3:length(Subfolders)
    Person = Subfolders(k).name;
    Files = dir(fullfile(InputDir + "/" + Person, "*.wav"));
    
    figure1 = figure('visible','off', 'Position', get(0, 'Screensize'));
    title(sprintf("DFT of voice of %s", Person));
    for j = 1:length(Files)
        [x, Fs] = audioread(InputDir + "/" + Person + "/" + Files(j).name);
        x = x';
        szx = size(x);
        if szx(1)>1
            x = (x(1, :)+x(2, :))/2;
        end
        len = min(len, length(x));
        Timing = Files(j).name;
        Timing = Timing(1:end-4);

        subplot(2, 3, j);
        plotDFT(x, Fs, length(x));
        xlim([-4e3 4e3]);
        title(sprintf("In the %s", Timing));

        fprintf(Person + " " + Files(j).name + " Processed\n")
    end
    
    FileDir = OutputDir + "/" + Person;
    if ~exist(FileDir, 'dir')
       mkdir(FileDir)
    end
    saveas(figure1, FileDir + "/" + "DFT", "png");
    
end

ET = [];
Column_Names = {};
Bands = {};
for i = 1:10
    Bands{end+1} = char(sprintf("B%d_%d", (i-1)*400, i*400));
end

All_fftx = [];

for k = 3:length(Subfolders)
    Person = Subfolders(k).name;
    Files = dir(fullfile(InputDir + "/" + Person, "*.wav"));
    
    figure1 = figure('visible','off', 'Position', get(0, 'Screensize'));
    title(sprintf("DFT of voice of %s", Person));
    for j = 1:length(Files)
        [x, Fs] = audioread(InputDir + "/" + Person + "/" + Files(j).name);
        x = x';
        szx = size(x);
        if szx(1)>1
            x = (x(1, :)+x(2, :))/2;
        end
        Timing = Files(j).name;
        Timing = Timing(1:end-4);
        Column_Names{end+1} = char(Person + "_" + Timing);
        
        Fx = fft(x);
        Fx = abs(Fx);
        All_fftx = [All_fftx; Fx(1:len)];
        
        All_y = [];
        All_ffty = [];
        for i = 1:10
            y = bandpass(x, [(i-1)*400+1 i*400], Fs);
            Fy = fft(y);
            Fy = abs(Fy);
            All_y = [All_y; y];
            All_ffty = [All_ffty; Fy(1:len)];
                        
            subplot(2, 5, i);
            plotDFT(y, Fs, length(x))
            xlim([-4e3 4e3]);
            title(sprintf("%d - %d", (i-1)*400, i*400));
        end
        
        FileDir = OutputDir + "/" + Person + "/" + "BandPass_Filtered";
        if ~exist(FileDir, 'dir')
           mkdir(FileDir)
        end
        saveas(figure1, FileDir + "/" + Timing, "png");
        
        ETP = [];
        for i=1:10
            ETP = [ETP bandpower(All_y(i, :))];
        end
        ET = [ET; ETP];
        
        Euc = zeros(10, 10);
        for a=1:10
            for b=1:10
                Euc(a, b) = sqrt(sum((All_ffty(a, :) - All_ffty(b, :)) .^ 2));
            end
        end
        EucTable = array2table(Euc, 'RowNames', Bands, 'VariableNames', Bands);
        writetable(EucTable, FileDir + "/" + Timing + "_Euclidean.csv");
        
        fprintf(Person + " " + Files(j).name + " Processed for Different Bandpass\n")
    end
    
end

EnergyTable = array2table(ET', 'VariableNames', Column_Names, 'RowNames', Bands);
writetable(EnergyTable, OutputDir + "/Energy_Values.csv");

Euc_Dist = zeros(length(Column_Names), length(Column_Names));
for i=1:length(Column_Names)
    for j=1:length(Column_Names)
        Euc_Dist(i, j) = sqrt(sum((All_fftx(i, :) - All_fftx(j, :)) .^ 2));
    end
end
EuclideanTable = array2table(Euc_Dist, 'VariableNames', Column_Names, 'RowNames', Column_Names);
writetable(EuclideanTable, OutputDir + "/Euclidean_Distance_Values.csv");

function plotDFT(x, Fs, N)

Fx = fft(x, N);
m = abs(Fx);
shifted_m = fftshift(m);

f = -Fs/2:Fs/N:Fs/2;
plot(f(1:N), 20*log10(shifted_m/N));
xlabel("Frequency (in Hz)"); ylabel("|X(f)/N|(dB)");

end