



sync      = bml_annot_read(['annot/' SUBJECT '_sync.txt']);
session   = bml_annot_read(['annot/' SUBJECT '_session.txt']);
electrode = bml_annot_read(['annot/' SUBJECT '_electrode.txt']);

NEW_PATH = '/Volumes/Nexus/DBS';
sync.folder = strrep(sync.folder, '\\136.142.16.9\Nexus\DBS', NEW_PATH);
sync.folder = strrep(sync.folder, '\', '/');

%loading neuroomega macro
cfg=[];
cfg.chantype = 'emg';
cfg.roi = sync(sync.filetype=="neuroomega.mat" & sync.chantype=="analog",:);
sync_no_macro=bml_sync_transfer_neuroomega_chantype(cfg);

cfg=[];
cfg.roi=sync_no_macro(end-1,:);
cfg.chantype = 'emg';
cfg.allow_missing = true;
cfg.electrode = electrode(electrode.type=="emg",:);
cfg.timetol = 0.001; %increase if necessary
no_emg = bml_load_continuous(cfg);

%%
cfg=[];
cfg.viewmode = 'vertical';
cfg.continuous = 'yes';
cfg.blocksize = 30;
ft_databrowser(cfg,no_emg)

%%%%%%%


electrode_emg_stim = electrode(electrode.type=="emg",:);

cfg=[];
cfg.roi=no_info(no_info.depth == -4.001,:)
cfg.chantype = 'emg';
cfg.filetype = 'neuroomega.mat';
cfg.allow_missing = true;
cfg.electrode = electrode_emg_stim;
cfg.electrode.starts = repmat(cfg.roi.starts(1),size(cfg.electrode.starts));
cfg.electrode.ends = repmat(cfg.roi.ends(end),size(cfg.electrode.starts));
cfg.timetol = 0.001; %increase if necessary
no_emg = bml_load_continuous(cfg);


%%
cfg=[];
cfg.viewmode = 'vertical';
cfg.continuous = 'yes';
cfg.blocksize = 180;
ft_databrowser(cfg,no_emg)
