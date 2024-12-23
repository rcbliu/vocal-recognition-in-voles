%%Fit sinFM function to each extracted vocal segment.
%%Prior to this, run USVSeg on the data along with our posthoc cleaning
%%code.
%%This then bandpass filters each of the vocalizations to only extract the
%%fundamenal contour (aka it excludes any harmonics) and fits the sinFM
%%function. 
%%I found that when harmonics were present, the fit would often jump
%%between harmonics, which was inherently problematic/incorrect.

clear
clc
close all

synology_drive = 'Z:';

plot_data = 1; %%plot the data (1) or not.
%%If you run the whole program with this on, you'll hate your life and it
%%will take 87 years (at minimum). This is only useful if you want to
%%ensure the accuracy of the fit.

base_path = fullfile(synology_drive, 'Megan\AdultRecs-FreeInteraction');

save_dir = fullfile(synology_drive, 'Megan\sinFM_GPU_Output\individVocs_fitSinFM\MakeNewFigs');
if ~exist(save_dir)
    mkdir(save_dir)
end

%%load in the list of files
[~,~,data] = xlsread(fullfile(synology_drive, 'Megan\AdultRecs-FreeInteraction\RecordingList'));
%%This data is organized such that each row corresponds to three separate
%%audio files (pretest, prepair, posttest). Read in all of that audio data.

if plot_data == 1
    fig1 = figure('color','w','position',[161.6667  198.3333  709.3333  396.0000]);
end

dur = -1; % ms. This was already here and I don't want to risk changing it. It just
%%gets overwritten later anyway

% optimization parameters
iterno = 2000; %%how many interations to try to optimize parameters.
d_F = .1;
d_Phi = 2*pi/100;

envcutoff = .001;
    
count = 1;
for dataset = 2:size(data,1) %%ignore header row
    tic
    
    for j = 1:3 %%for pretest (Day 0), prepair (Day 1), and postpair (Day 9):
        if j == 1
            %%first convert date
            date_tmp = data{dataset,6};
            label_tmp = '30MinMFDyads-preTest';
            audio_fname = data{dataset,7};
        elseif j == 2
            date_tmp = data{dataset,8};
            label_tmp = '30MinMFDyads-prePairing';
            audio_fname = data{dataset,9};
        elseif j == 3
            date_tmp = data{dataset,10};
            label_tmp = '30MinMFDyads-postPairing';
            audio_fname = data{dataset,11};
        end
        %%get the relevant info from the file path
        backslash = strfind(date_tmp, '/');
        month = date_tmp(1:backslash(1)-1);
        day = date_tmp(backslash(1)+1:backslash(2)-1);
        year = date_tmp(backslash(2)+1:end);
        %%This is annoying and convoluted but alas. Here we are.
        date_str = sprintf('%02d%02d%04d', str2double(month),str2double(day),str2double(year));
        base_path_tmp = fullfile(base_path, label_tmp, date_str, 'Audio');
        audio_path = base_path_tmp;
        contour_path = fullfile(base_path_tmp, 'ContourOutput');
        
        if strcmp(date_str, '05232021')%%I messed up this file. Ignore it for now.
            continue
        end

        %%if output already exists, skip
        if exist(fullfile(save_dir, sprintf('all_call_data_%s.mat', audio_fname)), 'file')
            continue
        end
        
        %%If file doesn't exist, run it.
        fprintf('Starting %s\n', audio_fname);
        try
            s = load(fullfile(contour_path, sprintf('SegmentData_%s', audio_fname)));
        catch
            s = load(fullfile(contour_path, 'Old', sprintf('SegmentData_%s', audio_fname)));
        end
        [audio_data,fs] = audioread(fullfile(audio_path, sprintf('%s.wav', audio_fname)));
        SegmentData = s.SegmentData; clear s
        
        %%preallocate space
        dur_list = NaN(length(SegmentData.StartTime),1);
        m_list = dur_list; %%slopes
        b_list = dur_list; %%offset freq
        A_list = dur_list; %%amplitude of freq modulation [Hz]
        F_list = dur_list; %%frequency of freq modulation (cyc/sec)
        Phi_list = dur_list; %%phi
        sum_diff = dur_list;  %%rmse
        
        for segment = 1:length(SegmentData.Contours)
            start_tmp = SegmentData.StartTime(segment)*fs;
            end_tmp = SegmentData.EndTime(segment)*fs;
            buffer = 0; %%add a buffer before/after each segment
            high_freq = nanmax(SegmentData.Contours{segment,1}(:,3));
            low_freq = nanmin(SegmentData.Contours{segment,1}(:,3));
            
            %%skip any vocs that fall below our freq threshold
            if high_freq < 22000
                continue
            end
            
            %%figure out bandpass freqs
            low_freq = low_freq-3000; %%add a frequency buffer
            high_freq = high_freq+3000;
            %%keep it within the confines of our freq ranges
            if low_freq < 22000
                low_freq = 22000;
            end
            if high_freq > 120000
                high_freq = 120000;
            end
            if low_freq > 120000
                continue
            end
            
            if isnan(SegmentData.StartTime)
                continue
            end
            
            audio_tmp = audio_data(start_tmp:end_tmp);
            %%suppress this warning
            w = warning('query','last');
            id = w.identifier;
            warning('off',id);
            
            filtered_audio = bandpass(audio_tmp, [low_freq high_freq], fs);
            
            if plot_data == 1
                figure(fig1);
                subplot(2,1,1);
                [s, f, t] = spectrogram(audio_tmp, 100, 1, 200, fs);
                surf(t, f, 20*log10(abs(s)), 'EdgeColor', 'none');
                axis xy;
                axis tight;
                colormap(jet); view(0,90);
                xlabel('Time (secs)');
                caxis([-40 max(audio_tmp)]);
                ylim([0 105000]);
                ytickvals = linspace(0, 105000, 4);
                set(gca, 'ytick', ytickvals);
                set(gca, 'yticklabel', ytickvals./1000);
                ylabel('Frequency (kHz)');
                xvals = xlim;
                set(gca, 'xtick', xvals);
                
                subplot(2,1,2);
                [s, f, t] = spectrogram(filtered_audio, 100, 1, 200, fs);
                surf(t, f, 20*log10(abs(s)), 'EdgeColor', 'none');
                axis xy;
                axis tight;
                colormap(jet); view(0,90);
                xlabel('Time (secs)');
                caxis([-40 max(filtered_audio)]);
                ylim([0 105000]);
                ytickvals = linspace(0, 105000, 4);
                set(gca, 'ytick', ytickvals);
                set(gca, 'yticklabel', ytickvals./1000);
                ylabel('Frequency (kHz)');
                xvals = xlim;
                set(gca, 'xtick', xvals);
            end
            sound = filtered_audio;
            clear filtered_audio audio_tmp
            
            clear smodel env fre s
            %%now fit the sinFM function to the newly filtered data
            [smodel,env,fre,s] = getcallhilb_tonemodel_direct(sound,'envcutoff', envcutoff, ...
                'samplerate',fs, ...
                'precallpts',0);
            
            if dur == -1
                calldur = max(find(env))/fs*1000;
            else
                calldur = dur;
            end
            
            %%The calldur code doesn't always work. Backup:
            if isempty(calldur)
                calldur = (SegmentData.EndTime(segment)-SegmentData.StartTime(segment))*1000;
            end
            
            n = find(env(1:round(fs * calldur / 1000)));
            
            Frn = fre(n)';
            
            % call is loaded now fit line to first '?' ms of frequency
            % vector using least squared linear regression
            Q = [n/fs, ones(length(n),1)];
            x = Q\Frn;
            
            m = x(1); % slope [Hz]/sec
            b = x(2); % onset
            Lin = m*n/fs + b;
            
            % compute residual
            rFrn = Frn - Lin;
            
            % fit sine wav to residual compute constants
            N = length(n);
            rFrs = sum(rFrn);
            
            % init parameters
            A = 1;
            DC = 0;
            F = 20; %Hz
            Phi = 0;
            
            % init vectors
            Sinn = sin(F*2*pi/fs*n-Phi);
            Cosn = cos(F*2*pi/fs*n-Phi);
            
            % init constants
            Sins = sum(Sinn);
            Coss = sum(Cosn);
            Sin2s = sum(Sinn.^2);
            SinrFrs = sum(Sinn.*rFrn);
            
            % %         cost = [];
            for i = 1:iterno
                % compute common vector for freq and phase partials
                xn = A*Sinn - rFrn + DC;
                
                % optimize frequency
                sp_F = sign(sum(n.*Cosn.*xn));
                F = F - sp_F*d_F*rand;
                
                % optimize phase
                sp_Phi = sign(-sum(Cosn.*xn));
                Phi = Phi - sp_Phi*d_Phi*rand;
                if Phi < 0
                    Phi = 2*pi + Phi; % keep positive phase
                end
                if Phi > 2*pi
                    Phi = Phi - 2*pi; % keep phase under 2pi
                end
                
                % update vectors
                Sinn = sin(F*2*pi/fs*n-Phi);
                Cosn = cos(F*2*pi/fs*n-Phi);
                
                % update constants
                Sins = sum(Sinn);
                Coss = sum(Cosn);
                Sin2s = sum(Sinn.^2);
                SinrFrs = sum(Sinn.*rFrn);
                
                A = SinrFrs / Sin2s; % don't take dc into account
            end
            
            if plot_data == 1
                fig2 = figure('color','w','position',[872.3333  199.0000  560.0000  394.0000]);
                subplot(2,2,1); hold off
                plot((1:length(smodel))/fs,smodel)
                hold on
                plot((1:length(env))/fs,env,'r')
                plot((1:length(env))/fs,env*-1,'r')
                hold off
                try
                    ylim([max(smodel)*2*-1, max(smodel*2)]);
                catch
                    continue
                end
                title(['Call-' num2str(segment)])
                ylabel('pressure')
                xlabel('time')
                
                subplot(2,2,2); hold off
                for i = find(env ~=0 )
                    plot(i/fs,fre(i),'x'); hold on
                end
                
                plot(n/fs,Lin + A*Sinn,'rx')
                hold off
                ylim([1e4, 10e4])
                xlim([0 n(end)/fs])
                title(['Call-' num2str(segment) ' spectrogram + fm fit'])
                ylabel('Frequency (kHz)');
                xlabel('Time (s)')
                
                subplot(2,2,3); hold off
                axis off
                
                text(.1,.9,['dur [ms]: ' num2str(calldur,3)])
                text(.1,.6,['A [Hz]: ' num2str(A,3)])
                text(.1,.5,['F [cyc/sec]: ' num2str(F,3)])
                text(.1,.7,['b [Hz]: ' num2str(b,3)])
                text(.1,.4,['Phi [deg]: ' num2str(Phi/(2*pi)*360,3)])
                text(.1,.8,['m [Hz/sec]: ' num2str(m,3)])
                
                subplot(2,2,4); hold off
                [s,f,t] = spectrogram(sound,100, 1, 200, fs);
                surf(t, f, 20*log10(abs(s)), 'EdgeColor', 'none');
                axis xy;
                axis tight;
                colormap(jet); view(0,90);
                xlabel('Time (s)');
                caxis([-40 max(audio_data)]);
                ylim([10000 100000]);
                yvals = ylim;
                ytickvals = linspace(yvals(1), yvals(2), 4);
                set(gca, 'ytick', ytickvals);
                set(gca, 'yticklabel', floor(ytickvals./1000));
                ylabel('Frequency (kHz)');
                xvals = xlim;
                xpos = (xvals(2)-xvals(1))*.7+xvals(1);
                xlim([0 n(end)/fs])
            end
            
            dur_list(segment,1) = calldur;
            m_list(segment,1) = m; %%slopes
            b_list(segment,1) = b; %%onset freq
            A_list(segment,1) = A; %%amplitude
            F_list(segment,1) = F; %%frequency modulation (Hz)
            Phi_list(segment,1) = Phi;
            
            if plot_data == 1
                pause(1.0);
                close(fig2)
            end
            
            %%calculate the octave difference RMSE from the contour 
            timestamps_tmp = SegmentData.Contours{segment,1}(:,1);
            timestamps_tmp = timestamps_tmp*1000; %%convert to ms
            timestamps_tmp = timestamps_tmp-timestamps_tmp(1);
            %%weird issue where the timestamps aren't in order. So this is
            %%awkward. Throw in a catch and reorganize if necessary
            try
                seg_freqs_interpolate = interp1(timestamps_tmp, SegmentData.Contours{segment,1}(:,3),...
                    [0:timestamps_tmp(end)/length(Lin):timestamps_tmp(end)], 'spline');
            catch
                data_tmp = cat(2,timestamps_tmp,SegmentData.Contours{segment,1}(:,3));
                %%weird issue with overlapping points. Delete them
                [B,I] = sort(data_tmp(:,1), 'ascend');
                data_tmp = data_tmp(I,:);
                diff_vals = diff(data_tmp(:,1));
                idx_del = find(diff_vals == 0);
                data_tmp(idx_del,:) = [];
                
                timestamps_tmp = data_tmp(:,1);
                frequencies_tmp = data_tmp(:,2);
                seg_freqs_interpolate = interp1(timestamps_tmp, frequencies_tmp,...
                    [0:timestamps_tmp(end)/length(Lin):timestamps_tmp(end)], 'spline');
                clear timestamps_tmp frequencies_tmp data_tmp
            end
            data_tmp = cat(2,seg_freqs_interpolate(1:end-1)',Lin+A*Sinn);
            
            %%now calculate the RMSE
            diff_val = sqrt(nanmean(log2(data_tmp(:,1)./data_tmp(:,2)).^2));
            sum_diff(segment,1) = diff_val;
            
            clear calldur m b A F Phi diff_val data_tmp timestamps_tmp seg_freqs_interpolate
        end
        
        % keep all Phi under 180 deg in param_mx
        Phi_list_deg = Phi_list/(2*pi)*360;
        A_list2 = A_list;
        for i  = 1:length(Phi_list_deg)
            if Phi_list_deg(i) >= 180
                Phi_list_deg(i) = Phi_list_deg(i) - 180; % keep phase under 180
                A_list2(i) = -A_list2(i);
            end
        end
        %%sum_diff = octave difference between USVSEG extracted contour and
        %%sinFM contour
        param_mx = [dur_list, A_list2, F_list, b_list, m_list, Phi_list_deg, sum_diff];
        
        save(fullfile(save_dir, sprintf('all_call_data_%s.mat', audio_fname)), ...
            'param_mx');

        fprintf('Done %i of %i files\n', (dataset-2)*3 + j, (size(data,1)-1)*3)
        toc
        clear sound
        clear dur_list m_list b_list A_list Phi_list
        clear F_list Phi_list_deg cost_list param_mx
    end
end