function random_real(){
    return Math.random();
}


function make_stateless_xrp(xrp_name, sampler, scorer, proposal_support){
    var sample = function(stats, hyperparams, args){
        var value = sampler.apply(null, args);
        return [value, stats, scorer(args, value)];
    };
    var incr_stats = function(value, stats, hyperparams, args){
        return [value, stats, scorer(args, value)];
    };
    var decr_stats = function(value, stats, hyperparams, args){
        return [value, stats, scorer(args, value)];
    };
    var proposer = null;
    if (proposal_support == []){
        proposer = null;
    }else{
        proposer = proposal_support[0];
    }

    var support = null;
    if (proposal_support.length < 2){
        support = null;
    }else{
        support = function(stats, hyperparams, args){return proposal_support[1](args);};
    }
    

    return church_make_xrp(xrp_name, sample, incr_stats, decr_stats, null, null, proposer, support);

}

var flip = make_stateless_xrp(
        "flip", 
        function(){
            return (arguments.length == 0)? (random_real() < 0.5): (random_real() < arguments[0]);
        },
        function(args, val){
            if (args.length == 0){
                return -Math.log(2.0);
            }else{
                return val? Math.log(args[0]): Math.log(1 - args[0]);
            }
        },
        [
            function(ops, old_val){
                var new_val = !old_val;
                return [new_val, 0.0, 0.0];
            },
            function(args){
                return [true, false];
            }
        ]
        );

