function tpcr_decomposition_script(start,send,step,subfolder)
    tframe=start:step:min(start+999,send);
    if (start-send)~=0
        tpcr_decomposition(tframe,subfolder);
    end
end