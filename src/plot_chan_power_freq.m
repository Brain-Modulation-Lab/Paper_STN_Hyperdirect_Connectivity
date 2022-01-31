function plot_chan_power_freq(cfg)


% frequency analysis
%
% plots power vs freq for x channels, all times
%
% ajorge 2018






%% --------------------------------------------------------------------------
% plot power vs freq for a selected channel



if cfg.flag == 2
    
    
    
    %--------------------------------------------------------------------------
    %matlab's pwelch method
    
    data = cfg.data;
    
    figure(330);
    [getlength] = pwelch(raw.trial{1}(1,1:stim_marker2(1)-1),[],[],[],1000);
    
    for ch = 3
        
        
        
        if max(raw.trial{1}(ch,1:stim_marker2(1)-1))<1500
            [pre_pxx(ch,:),pre_f(ch,:)] = pwelch(raw.trial{1}(ch,1:stim_marker2(1)-1),[],[],[],1000);
            subplot(2,3,1)
            semilogx(pre_f(ch,:),smooth(10*log10(pre_pxx(ch,:))),'Color', [200 200 200]./256)
            hold on

        else
            pre_pxx(ch,1:getlength) = NaN;
            pre_f(ch,1:getlength)   = NaN;
        end 
    end

    %prestim mean
    clear pre_pxx_mean
    title('Raw: Pre stim - all ch')
    pre_pxx_mean = nanmedian(pre_pxx(:,:),1);
    subplot(2,3,1)
    semilogx(pre_f(1,:),smooth(10*log10(pre_pxx_mean)),'k','LineWidth',2)
    xlim([10^0 10^2.699])
    ylim([-40 70])
    %set(gca,'xticklabel',{[]}) 
    hold on


   
    %{
    %--------------------------------------------------------------------------
    %fieldtrip method
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.output      = 'pow';
    cfg.channel     = 'lfp 3';
    cfg.trials      = 'all';
    cfg.taper       = 'hanning';
    cfg.foi         = 2:0.01:80;                           % analysis 2 to 30 Hz in steps of 2 Hz 
    cfg.t_ftimwin   = ones(length(cfg.foi),1).*0.05;        % length of time window = 0.5 sec
    %cfg.t_ftimwin  = 7./cfg.foi;  % 7 cycles per time window
    cfg.toi         = 'all';                       % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
    %cfg.pad        = 'nextpow2';

    hann = ft_freqanalysis(cfg, raw);
    figure(101)
    semilogx(hann.freq, smooth(hann.powspctrm))
    %}





end





%% --------------------------------------------------------------------------
% plot power vs freq (pwelch) on all channels on all times
if cfg.flag == 1

    
    
    %--------------------------------------------------------------------------
    % Frequency analysis Raw - all channels

    clear pre_pxx pre_f

    %prestim
    figure(330);
    [getlength] = pwelch(raw.trial{1}(1,1:stim_marker2(1)-1),[],[],[],1000);
    for ch = 1:63
        if max(raw.trial{1}(ch,1:stim_marker2(1)-1))<1500
            [pre_pxx(ch,:),pre_f(ch,:)] = pwelch(raw.trial{1}(ch,1:stim_marker2(1)-1),[],[],[],1000);
            subplot(2,3,1)
            semilogx(pre_f(ch,:),smooth(10*log10(pre_pxx(ch,:))),'Color', [200 200 200]./256)
            hold on

        else
            pre_pxx(ch,1:getlength) = NaN;
            pre_f(ch,1:getlength)   = NaN;
        end
    end

    %prestim mean
    clear pre_pxx_mean
    title('Raw: Pre stim - all ch')
    pre_pxx_mean = nanmedian(pre_pxx(:,:),1);
    subplot(2,3,1)
    semilogx(pre_f(1,:),smooth(10*log10(pre_pxx_mean)),'k','LineWidth',2)
    xlim([10^0 10^2.699])
    ylim([-40 70])
    %set(gca,'xticklabel',{[]}) 
    hold on

    % find data for all stim sets
    clear stim_pxx stim_f
    for set = 1:stim_marker2(end,2)
        indx = find(stim_marker2(:,2)==set);
        begg = stim_marker2(indx(1),1);
        endd = stim_marker2(indx(end),1);
        for ch = 1:63           
                [stim_pxx(ch,:),stim_f(ch,:)] = pwelch(raw.trial{1}(ch,begg:endd),[],[],[],1000);
                subplot(2,3,4)
                semilogx(stim_f(ch,:),smooth(10*log10(stim_pxx(ch,:))),'Color', [200 200 200]./256)
                hold on
        end

    end


    %stim mean
    title('Raw: During stim - all ch')
    clear stim_pxx_mean
    stim_pxx_mean = nanmedian(stim_pxx(:,:),1);
    subplot(2,3,4)
    semilogx(stim_f(1,:),smooth(10*log10(stim_pxx_mean)),'k','LineWidth',2)
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB/Hz)')
    xlim([10^0 10^2.699])
    ylim([-40 70])

    hold on










    %--------------------------------------------------------------------------
    %{
    %% Frequency analysis - only if ch < iq_voltage75 && ch<iq_voltage25

    iq_upper = 0.75;
    iq_lower = 0.25;

    %separate into outlier channels and the rest
    median_voltagepre   = nanmedian(pre_voltage_mean(:,pt),1);
    iq_voltage25pre    = quantile(pre_voltage_mean(:,pt),iq_lower);
    iq_voltage75pre     = quantile(pre_voltage_mean(:,pt),iq_upper);

    median_voltagestim   = nanmedian(set_voltage_mean2(:,pt),1);
    iq_voltage25stim     = quantile(set_voltage_mean2(:,pt),iq_lower);
    iq_voltage75stim     = quantile(set_voltage_mean2(:,pt),iq_upper);


    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % IQ VOLTAGE < 75 QUANTILE 

    %prestim
    figure(330);

    clear pre_pxx pre_f
    [getlength] = length(pwelch(raw.trial{1}(1,1:stim_marker2(1)-1),[],[],[],1000));

    for ch = 1:63
       % if max(raw.trial{1}(ch,1:stim_marker2(1)-1))<1000
        if pre_voltage_mean(ch,pt) < iq_voltage75pre
            [pre_pxx(ch,:,pt),pre_f(ch,:,pt)] = pwelch(raw.trial{1}(ch,1:stim_marker2(1)-1),[],[],[],1000);
            subplot(2,3,2)
            semilogx(pre_f(ch,:,pt),smooth(10*log10(pre_pxx(ch,:,pt))),'Color', [200 200 200]./256)
            hold on

        else
            lengthh = length(raw.trial{1}(ch,1:stim_marker2(1)-1));
            pre_pxx(ch,1:getlength,pt) = NaN;
            pre_f(ch,1:getlength,pt)   = NaN;
        end
    end

    %prestim mean
    clear pre_pxx_mean
    title('Raw: Pre stim - ch>75q')
    pre_pxx_mean = nanmedian(pre_pxx(:,:,pt),1);
    subplot(2,3,2)
    semilogx(pre_f(1,:,pt),smooth(10*log10(pre_pxx_mean)),'k','LineWidth',2)
    xlim([10^0 10^2.699])   %10^2.699 = 500 Hz
    ylim([-40 70])
    hold on



    % find data for all stim sets
    for set = 1:stim_marker2(end,2)
        indx = find(stim_marker2(:,2)==set);
        begg = stim_marker2(indx(1),1);
        endd = stim_marker2(indx(end),1);
        for ch = 1:63           
            if set_voltage_mean2(ch,pt) < iq_voltage75stim
                [stim_pxx(ch,:,pt),stim_f(ch,:,pt)] = pwelch(raw.trial{1}(ch,begg:endd),[],[],[],1000);
                subplot(2,3,5)
                semilogx(stim_f(ch,:,pt),smooth(10*log10(stim_pxx(ch,:,pt))),'Color', [200 200 200]./256)
                hold on

            else
                lengthh = length(raw.trial{1}(ch,1:stim_marker2(1)-1));
                stim_pxx(ch,1:getlength,pt) = NaN;
                stim_f(ch,1:getlength,pt)   = NaN;
            end


        end

    end


    %stim mean
    clear stim_pxx_mean
    title('Raw: During stim - ch>75q')
    stim_pxx_mean = nanmedian(stim_pxx(:,:,pt),1);
    subplot(2,3,5)
    semilogx(stim_f(1,:,pt),smooth(10*log10(stim_pxx_mean)),'k','LineWidth',2)
    xlim([10^0 10^2.699])
        ylim([-40 70])
    hold on



    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % IQ VOLTAGE < 25 QUANTILE 

    clear pre_pxx pre_f
    %pre_pxx(1:63,1:getlength,pt) = NaN*ones(1:63,1:getlength);

    getlength = length(pwelch(raw.trial{1}(ch,1:stim_marker2(1)-1),[],[],[],1000));
    for ch = 1:63
       % if max(raw.trial{1}(ch,1:stim_marker2(1)-1))<1000
        if pre_voltage_mean(ch,pt) < iq_voltage25pre
            [pre_pxx(ch,:,pt),pre_f(ch,:,pt)] = pwelch(raw.trial{1}(ch,1:stim_marker2(1)-1),[],[],[],1000);
            subplot(2,3,3)
            semilogx(pre_f(ch,:,pt),smooth(10*log10(pre_pxx(ch,:,pt))),'Color', [200 200 200]./256)
            hold on
            chforplot=ch;

        else
            pre_pxx(ch,1:getlength,pt) = NaN;
            pre_f(ch,1:getlength,pt)   = NaN;

        end
    end

    %prestim mean
    clear pre_pxx_mean
    title('Raw: Pre stim - ch>25q')
    pre_pxx_mean = nanmedian(10*log10(pre_pxx(:,:,pt)));
    subplot(2,3,3)
    semilogx(pre_f(chforplot,:,pt),smooth((pre_pxx_mean)),'k','LineWidth',2)
    xlim([10^0 10^2.699])   %10^2.699 = 500 Hz
    ylim([-40 70])
    hold on

    clear stim_pxx stim_f
    % find data for all stim sets
    counter3=1;
    for set = 1:stim_marker2(end,2)
        indx = find(stim_marker2(:,2)==set);
        begg = stim_marker2(indx(1),1);
        endd = stim_marker2(indx(end),1);
        for ch = 1:63           
            if set_voltage_mean2(ch,pt) < iq_voltage25stim
                counter3=counter3+1;
                [stim_pxx(ch,:,pt),stim_f(ch,:,pt)] = pwelch(raw.trial{1}(ch,begg:endd),[],[],[],1000);
                subplot(2,3,6)
                semilogx(stim_f(ch,:,pt),smooth(10*log10(stim_pxx(ch,:,pt))),'Color', [200 200 200]./256)
                hold on

            else
                lengthh = length(raw.trial{1}(ch,1:stim_marker2(1)-1));
                stim_pxx(ch,1:getlength,pt) = NaN;
                stim_f(ch,1:getlength,pt)   = NaN;
            end


        end

    end


    %stim mean
    clear stim_pxx_mean
    title('Raw: During stim - ch>25q')
    stim_pxx_mean = nanmedian(stim_pxx(:,:,pt),1);
    subplot(2,3,6)
    semilogx(stim_f(1,:,pt),smooth(10*log10(stim_pxx_mean)),'k','LineWidth',2)
    xlim([10^0 10^2.699])
        ylim([-40 70])
    hold on

    %}

end
    