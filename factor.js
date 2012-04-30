function Factor(address, args, value, factor_function, ticks){
    this.address = address;
    this.args = args;
    this.value = value;
    this.factor_function = factor_function;
    this.ticks = ticks;
}

function make_factor(factor_function){
    return function(){
        var args = arguments;

        var new_val = 0;
        var factor_address = getCurrentAddress();

        var update_fx = function(factor_instance){
            new_val = factor_function.apply(null, args);
            var last_tick = false;
            if (found(factor_instance)){
                last_tick = factor_instance.ticks[0];
            }

            var new_factor = new Factor(factor_address, args, new_val, factor_function, [getCurrentStore().tick, last_tick]);
            getCurrentStore().score = getCurrentStore().score + new_val;
            return new_factor;
        };

        update_addbox(getCurrentStore().factors, factor_address, update_fx);
        return new_val;
    }

}


