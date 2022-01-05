
function stimulus =  load_sf_stimuli()

stim_run_fns ={'con', 'exp','clw','ccw'};
 
for i=1:length(stim_run_fns)
    stim_fn = sprintf('resources/stimuli/retmap%s.mat',stim_run_fns{i});
    stimstr= load(stim_fn);
    stimulus{i} = stimstr.stimuli; 
    
  
end   
 
end