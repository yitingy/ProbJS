function mh_query(samples, lag, normal_form_proc){
    curr_store_idx = 0;
    var store = getCurrentStore(); 
    var count = 0;
    var result_samples = [];

    proposal = function(store){
        var xrp_address_keys = Object.keys(store.xrp_draws);  
        if (xrp_address_keys.length == 0){
            return 0;
        }
        var selected_idx = Math.floor(Math.random()*xrp_address_keys.length);
        var new_val_fw_bw = store.xrp_draws[xrp_address_keys[selected_idx]].proposer_thunk();
        var new_val = new_val_fw_bw[0];
        store.xrp_draws[xrp_address_keys[selected_idx]].value = new_val;
        return (new_val_fw_bw[2] - new_val_fw_bw[1]);
    };

    while(count < samples){

        var old_score = getCurrentStore().score;
        var proposal_store = getCurrentStore().clone();
        curr_store_idx = 1;
        stores[curr_store_idx] = proposal_store;
        proposal_store.score = 0;
        proposal_store.tick += 1;
        var proposal_correction = proposal(proposal_store);
        curr_sample = normal_form_proc();

        console.log("proposed_sample");
        console.log(curr_sample);
        console.log("old_score");
        console.log(old_score);
        console.log("new_score");
        console.log(proposal_store.score);
        var accept = Math.min(0, proposal_store.score - old_score + proposal_correction);
        if (Math.log(random_real()) < accept || count == 0){
            console.log("accept ===");
            //console.log(curr_sample);
            //console.log("proposal_store:");
            //console.log(proposal_store);
            //console.log("store:");
            //console.log(store);
            stores[0] = proposal_store;
            result_samples.push(curr_sample);
        }else{
            stores[0].tick += 1;
            result_samples.push(result_samples[result_samples.length-1]);
        }
        curr_store_idx = 0;
        count += 1;

    }
    return result_samples;
}

