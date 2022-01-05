
function stimulus =  load_hcp_stimuli()


stim_run_fns ={'RETCCW', 'RETCW','RETEXP','RETCON', 'RETBAR1', 'RETBAR2'};
for i=1:6
    stim_fn = sprintf('pRF-decoder/apertures/%ssmall.mat',stim_run_fns{i});
    stimstr= load(stim_fn);
    stimulus{i} = stimstr.stim;
end 

end