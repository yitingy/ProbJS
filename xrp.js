Object.prototype.clone = function() {
  var newObj = (this instanceof Array) ? [] : {};
  for (i in this) {
    if (i == 'clone') continue;
    if (this[i] && typeof this[i] == "object") {
      newObj[i] = this[i].clone();
    } else newObj[i] = this[i]
  } return newObj;
};


var default_store = {
        tick: 0,
        xrp_draws: {},
        xrp_stats: {},
        factors:{},
        score: 0
    };

var stores = [default_store, default_store];
var curr_store_idx = 0;


function getCurrentStore(){
    return stores[curr_store_idx];
}

function getCurrentAddress(){
    var trace = printStackTrace();
    var clean_trace = [];
    for (var i = 5; i < trace.length-5; i++){
        clean_trace.push(trace[i]);
    }
    return clean_trace;
}

function XRP_draw(address, value, xrp_name, proposer_thunk, ticks, score, support){
    this.address = address;
    this.value = value;
    this.xrp_name = xrp_name;
    this.proposer_thunk = proposer_thunk;
    this.ticks = ticks;
    this.score = score;
    this.support = support;
}

function update_addbox(stats, address, update_fx){
    stats[address] = update_fx(stats[address]);
}

function not_found(prop){
    if (typeof(prop) === "undefined"){
        return true;
    } else{
        return false;
    }
}

function found(prop){
    return !not_found(prop);
}

function church_make_xrp(xrp_name, sample, incr_stats, decr_stats, init_stats, hyperparams, _proposer, support){
    var xrp_address = getCurrentAddress();
    var update_fx = function(stats){
        if (not_found(stats) || stats.tick != getCurrentStore().tick){
            return [init_stats, getCurrentStore().tick];
        }else{
            return stats;
        }
    };
    update_addbox(getCurrentStore().xrp_stats, xrp_address, update_fx);

    var proposer = _proposer;
    if (proposer == null){
        proposer = function(operands, old_value){
            var dec = decr_stats(old_value, getCurrentStore().xrp_stats[xrp_address][0], hyperparams, operands);
            var decstats = dec[1];
            var decscore = dec[2];
            var inc = sample(decstats, hyperparams, operands);
            var proposal_value = inc[0];
            var incscore = inc[2];
            return [proposal_value, incscore, decscore];
        };
    }

   
    return function(){
        var args = arguments;
        var xrp_draw_address = getCurrentAddress();
        var new_val = 0;
        var update_xrp_fx = function(xrp_draw){
            if (found(xrp_draw) && (getCurrentStore().tick == xrp_draw.ticks[0])) {
                new_val = xrp_draw.value;
                getCurrentStore().score = getCurrentStore().score + xrp_draw.score;
                return xrp_draw;
            }else{
                var stats = getCurrentStore().xrp_stats[xrp_address][0];
                var support_vals = [];
                if(support != null){
                    support_vals = support(stats, hyperparams, args);
                }
                var tmp = null; 
                if(not_found(xrp_draw)){
                    tmp = sample(stats, hyperparams, args);
                }else{
                    tmp = incr_stats(xrp_draw.value, stats, hyperparams, args);
                }
                var value = tmp[0];
                var new_stats = [tmp[1], getCurrentStore().tick];
                var incr_score = tmp[2];
                var last_tick = false;
                if (found(xrp_draw)){
                    last_tick = xrp_draw.ticks[0];
                }
                var local_proposer_thunk = function(state){
                    var returned_val = proposer(args, value);
                    return returned_val;
                };
                var new_xrp_draw = new XRP_draw(xrp_draw_address, value, xrp_name, local_proposer_thunk, [getCurrentStore().tick, last_tick], incr_score, support_vals);
                new_val = value;
                getCurrentStore().xrp_stats[xrp_address] = new_stats;
                getCurrentStore().score = getCurrentStore().score + incr_score;
                return new_xrp_draw;
            }

        };
        update_addbox(getCurrentStore().xrp_draws, xrp_draw_address, update_xrp_fx);
        return new_val;
    };

}
