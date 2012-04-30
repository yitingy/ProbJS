// Domain Public by Eric Wendelin http://eriwen.com/ (2008)
//                  Luke Smith http://lucassmith.name/ (2008)
//                  Loic Dachary <loic@dachary.org> (2008)
//                  Johan Euphrosine <proppy@aminche.com> (2008)
//                  Oyvind Sean Kinsey http://kinsey.no/blog (2010)
//                  Victor Homyakov <victor-homyakov@users.sourceforge.net> (2010)

/**
 * Main function giving a function stack trace with a forced or passed in Error
 *
 * @cfg {Error} e The error to create a stacktrace from (optional)
 * @cfg {Boolean} guess If we should try to resolve the names of anonymous functions
 * @return {Array} of Strings with functions, lines, files, and arguments where possible
 */
Error.stackTraceLimit = 1000;

function printStackTrace(options) {
    options = options || {guess: true};
    var ex = options.e || null, guess = !!options.guess;
    var p = new printStackTrace.implementation(), result = p.run(ex);
    return (guess) ? p.guessAnonymousFunctions(result) : result;
}

printStackTrace.implementation = function() {
};

printStackTrace.implementation.prototype = {
    /**
     * @param {Error} ex The error to create a stacktrace from (optional)
     * @param {String} mode Forced mode (optional, mostly for unit tests)
     */
    run: function(ex, mode) {
        ex = ex || this.createException();
        // examine exception properties w/o debugger
        //for (var prop in ex) {alert("Ex['" + prop + "']=" + ex[prop]);}
        mode = mode || this.mode(ex);
        if (mode === 'other') {
            return this.other(arguments.callee);
        } else {
            return this[mode](ex);
        }
    },

    createException: function() {
        try {
            this.undef();
        } catch (e) {
            return e;
        }
    },

    /**
     * Mode could differ for different exception, e.g.
     * exceptions in Chrome may or may not have arguments or stack.
     *
     * @return {String} mode of operation for the exception
     */
    mode: function(e) {
        if (e['arguments'] && e.stack) {
            return 'chrome';
        } else if (typeof e.message === 'string' && typeof window !== 'undefined' && window.opera) {
            // e.message.indexOf("Backtrace:") > -1 -> opera
            // !e.stacktrace -> opera
            if (!e.stacktrace) {
                return 'opera9'; // use e.message
            }
            // 'opera#sourceloc' in e -> opera9, opera10a
            if (e.message.indexOf('\n') > -1 && e.message.split('\n').length > e.stacktrace.split('\n').length) {
                return 'opera9'; // use e.message
            }
            // e.stacktrace && !e.stack -> opera10a
            if (!e.stack) {
                return 'opera10a'; // use e.stacktrace
            }
            // e.stacktrace && e.stack -> opera10b
            if (e.stacktrace.indexOf("called from line") < 0) {
                return 'opera10b'; // use e.stacktrace, format differs from 'opera10a'
            }
            // e.stacktrace && e.stack -> opera11
            return 'opera11'; // use e.stacktrace, format differs from 'opera10a', 'opera10b'
        } else if (e.stack) {
            return 'firefox';
        }
        return 'other';
    },

    /**
     * Given a context, function name, and callback function, overwrite it so that it calls
     * printStackTrace() first with a callback and then runs the rest of the body.
     *
     * @param {Object} context of execution (e.g. window)
     * @param {String} functionName to instrument
     * @param {Function} function to call with a stack trace on invocation
     */
    instrumentFunction: function(context, functionName, callback) {
        context = context || window;
        var original = context[functionName];
        context[functionName] = function instrumented() {
            callback.call(this, printStackTrace().slice(4));
            return context[functionName]._instrumented.apply(this, arguments);
        };
        context[functionName]._instrumented = original;
    },

    /**
     * Given a context and function name of a function that has been
     * instrumented, revert the function to it's original (non-instrumented)
     * state.
     *
     * @param {Object} context of execution (e.g. window)
     * @param {String} functionName to de-instrument
     */
    deinstrumentFunction: function(context, functionName) {
        if (context[functionName].constructor === Function &&
                context[functionName]._instrumented &&
                context[functionName]._instrumented.constructor === Function) {
            context[functionName] = context[functionName]._instrumented;
        }
    },

    /**
     * Given an Error object, return a formatted Array based on Chrome's stack string.
     *
     * @param e - Error object to inspect
     * @return Array<String> of function calls, files and line numbers
     */
    chrome: function(e) {
        var stack = (e.stack + '\n').replace(/^\S[^\(]+?[\n$]/gm, '').
          replace(/^\s+(at eval )?at\s+/gm, '').
          replace(/^([^\(]+?)([\n$])/gm, '{anonymous}()@$1$2').
          replace(/^Object.<anonymous>\s*\(([^\)]+)\)/gm, '{anonymous}()@$1').split('\n');
        stack.pop();
        return stack;
    },

    /**
     * Given an Error object, return a formatted Array based on Firefox's stack string.
     *
     * @param e - Error object to inspect
     * @return Array<String> of function calls, files and line numbers
     */
    firefox: function(e) {
        return e.stack.replace(/(?:\n@:0)?\s+$/m, '').replace(/^\(/gm, '{anonymous}(').split('\n');
    },

    opera11: function(e) {
        // "Error thrown at line 42, column 12 in <anonymous function>() in file://localhost/G:/js/stacktrace.js:\n"
        // "Error thrown at line 42, column 12 in <anonymous function: createException>() in file://localhost/G:/js/stacktrace.js:\n"
        // "called from line 7, column 4 in bar(n) in file://localhost/G:/js/test/functional/testcase1.html:\n"
        // "called from line 15, column 3 in file://localhost/G:/js/test/functional/testcase1.html:\n"
        var ANON = '{anonymous}', lineRE = /^.*line (\d+), column (\d+)(?: in (.+))? in (\S+):$/;
        var lines = e.stacktrace.split('\n'), result = [];

        for (var i = 0, len = lines.length; i < len; i += 2) {
            var match = lineRE.exec(lines[i]);
            if (match) {
                var location = match[4] + ':' + match[1] + ':' + match[2];
                var fnName = match[3] || "global code";
                fnName = fnName.replace(/<anonymous function: (\S+)>/, "$1").replace(/<anonymous function>/, ANON);
                result.push(fnName + '@' + location + ' -- ' + lines[i + 1].replace(/^\s+/, ''));
            }
        }

        return result;
    },

    opera10b: function(e) {
        // "<anonymous function: run>([arguments not available])@file://localhost/G:/js/stacktrace.js:27\n" +
        // "printStackTrace([arguments not available])@file://localhost/G:/js/stacktrace.js:18\n" +
        // "@file://localhost/G:/js/test/functional/testcase1.html:15"
        var lineRE = /^(.*)@(.+):(\d+)$/;
        var lines = e.stacktrace.split('\n'), result = [];

        for (var i = 0, len = lines.length; i < len; i++) {
            var match = lineRE.exec(lines[i]);
            if (match) {
                var fnName = match[1]? (match[1] + '()') : "global code";
                result.push(fnName + '@' + match[2] + ':' + match[3]);
            }
        }

        return result;
    },

    /**
     * Given an Error object, return a formatted Array based on Opera 10's stacktrace string.
     *
     * @param e - Error object to inspect
     * @return Array<String> of function calls, files and line numbers
     */
    opera10a: function(e) {
        // "  Line 27 of linked script file://localhost/G:/js/stacktrace.js\n"
        // "  Line 11 of inline#1 script in file://localhost/G:/js/test/functional/testcase1.html: In function foo\n"
        var ANON = '{anonymous}', lineRE = /Line (\d+).*script (?:in )?(\S+)(?:: In function (\S+))?$/i;
        var lines = e.stacktrace.split('\n'), result = [];

        for (var i = 0, len = lines.length; i < len; i += 2) {
            var match = lineRE.exec(lines[i]);
            if (match) {
                var fnName = match[3] || ANON;
                result.push(fnName + '()@' + match[2] + ':' + match[1] + ' -- ' + lines[i + 1].replace(/^\s+/, ''));
            }
        }

        return result;
    },

    // Opera 7.x-9.2x only!
    opera9: function(e) {
        // "  Line 43 of linked script file://localhost/G:/js/stacktrace.js\n"
        // "  Line 7 of inline#1 script in file://localhost/G:/js/test/functional/testcase1.html\n"
        var ANON = '{anonymous}', lineRE = /Line (\d+).*script (?:in )?(\S+)/i;
        var lines = e.message.split('\n'), result = [];

        for (var i = 2, len = lines.length; i < len; i += 2) {
            var match = lineRE.exec(lines[i]);
            if (match) {
                result.push(ANON + '()@' + match[2] + ':' + match[1] + ' -- ' + lines[i + 1].replace(/^\s+/, ''));
            }
        }

        return result;
    },

    // Safari, IE, and others
    other: function(curr) {
        var ANON = '{anonymous}', fnRE = /function\s*([\w\-$]+)?\s*\(/i, stack = [], fn, args, maxStackSize = 10;
        while (curr && curr['arguments'] && stack.length < maxStackSize) {
            fn = fnRE.test(curr.toString()) ? RegExp.$1 || ANON : ANON;
            args = Array.prototype.slice.call(curr['arguments'] || []);
            stack[stack.length] = fn + '(' + this.stringifyArguments(args) + ')';
            curr = curr.caller;
        }
        return stack;
    },

    /**
     * Given arguments array as a String, subsituting type names for non-string types.
     *
     * @param {Arguments} object
     * @return {Array} of Strings with stringified arguments
     */
    stringifyArguments: function(args) {
        var result = [];
        var slice = Array.prototype.slice;
        for (var i = 0; i < args.length; ++i) {
            var arg = args[i];
            if (arg === undefined) {
                result[i] = 'undefined';
            } else if (arg === null) {
                result[i] = 'null';
            } else if (arg.constructor) {
                if (arg.constructor === Array) {
                    if (arg.length < 3) {
                        result[i] = '[' + this.stringifyArguments(arg) + ']';
                    } else {
                        result[i] = '[' + this.stringifyArguments(slice.call(arg, 0, 1)) + '...' + this.stringifyArguments(slice.call(arg, -1)) + ']';
                    }
                } else if (arg.constructor === Object) {
                    result[i] = '#object';
                } else if (arg.constructor === Function) {
                    result[i] = '#function';
                } else if (arg.constructor === String) {
                    result[i] = '"' + arg + '"';
                } else if (arg.constructor === Number) {
                    result[i] = arg;
                }
            }
        }
        return result.join(',');
    },

    sourceCache: {},

    /**
     * @return the text from a given URL
     */
    ajax: function(url) {
        var req = this.createXMLHTTPObject();
        if (req) {
            try {
                req.open('GET', url, false);
                //req.overrideMimeType('text/plain');
                //req.overrideMimeType('text/javascript');
                req.send(null);
                //return req.status == 200 ? req.responseText : '';
                return req.responseText;
            } catch (e) {
            }
        }
        return '';
    },

    /**
     * Try XHR methods in order and store XHR factory.
     *
     * @return <Function> XHR function or equivalent
     */
    createXMLHTTPObject: function() {
        var xmlhttp, XMLHttpFactories = [
            function() {
                return new XMLHttpRequest();
            }, function() {
                return new ActiveXObject('Msxml2.XMLHTTP');
            }, function() {
                return new ActiveXObject('Msxml3.XMLHTTP');
            }, function() {
                return new ActiveXObject('Microsoft.XMLHTTP');
            }
        ];
        for (var i = 0; i < XMLHttpFactories.length; i++) {
            try {
                xmlhttp = XMLHttpFactories[i]();
                // Use memoization to cache the factory
                this.createXMLHTTPObject = XMLHttpFactories[i];
                return xmlhttp;
            } catch (e) {
            }
        }
    },

    /**
     * Given a URL, check if it is in the same domain (so we can get the source
     * via Ajax).
     *
     * @param url <String> source url
     * @return False if we need a cross-domain request
     */
    isSameDomain: function(url) {
        return typeof location !== "undefined" && url.indexOf(location.hostname) !== -1; // location may not be defined, e.g. when running from nodejs.
    },

    /**
     * Get source code from given URL if in the same domain.
     *
     * @param url <String> JS source URL
     * @return <Array> Array of source code lines
     */
    getSource: function(url) {
        // TODO reuse source from script tags?
        if (!(url in this.sourceCache)) {
            this.sourceCache[url] = this.ajax(url).split('\n');
        }
        return this.sourceCache[url];
    },

    guessAnonymousFunctions: function(stack) {
        for (var i = 0; i < stack.length; ++i) {
            var reStack = /\{anonymous\}\(.*\)@(.*)/,
                reRef = /^(.*?)(?::(\d+))(?::(\d+))?(?: -- .+)?$/,
                frame = stack[i], ref = reStack.exec(frame);

            if (ref) {
                var m = reRef.exec(ref[1]);
                if (m) { // If falsey, we did not get any file/line information
                    var file = m[1], lineno = m[2], charno = m[3] || 0;
                    if (file && this.isSameDomain(file) && lineno) {
                        var functionName = this.guessAnonymousFunction(file, lineno, charno);
                        stack[i] = frame.replace('{anonymous}', functionName);
                    }
                }
            }
        }
        return stack;
    },

    guessAnonymousFunction: function(url, lineNo, charNo) {
        var ret;
        try {
            ret = this.findFunctionName(this.getSource(url), lineNo);
        } catch (e) {
            ret = 'getSource failed with url: ' + url + ', exception: ' + e.toString();
        }
        return ret;
    },

    findFunctionName: function(source, lineNo) {
        // FIXME findFunctionName fails for compressed source
        // (more than one function on the same line)
        // TODO use captured args
        // function {name}({args}) m[1]=name m[2]=args
        var reFunctionDeclaration = /function\s+([^(]*?)\s*\(([^)]*)\)/;
        // {name} = function ({args}) TODO args capture
        // /['"]?([0-9A-Za-z_]+)['"]?\s*[:=]\s*function(?:[^(]*)/
        var reFunctionExpression = /['"]?([0-9A-Za-z_]+)['"]?\s*[:=]\s*function\b/;
        // {name} = eval()
        var reFunctionEvaluation = /['"]?([0-9A-Za-z_]+)['"]?\s*[:=]\s*(?:eval|new Function)\b/;
        // Walk backwards in the source lines until we find
        // the line which matches one of the patterns above
        var code = "", line, maxLines = Math.min(lineNo, 20), m, commentPos;
        for (var i = 0; i < maxLines; ++i) {
            // lineNo is 1-based, source[] is 0-based
            line = source[lineNo - i - 1];
            commentPos = line.indexOf('//');
            if (commentPos >= 0) {
                line = line.substr(0, commentPos);
            }
            // TODO check other types of comments? Commented code may lead to false positive
            if (line) {
                code = line + code;
                m = reFunctionExpression.exec(code);
                if (m && m[1]) {
                    return m[1];
                }
                m = reFunctionDeclaration.exec(code);
                if (m && m[1]) {
                    //return m[1] + "(" + (m[2] || "") + ")";
                    return m[1];
                }
                m = reFunctionEvaluation.exec(code);
                if (m && m[1]) {
                    return m[1];
                }
            }
        }
        return '(?)';
    }
};
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

function repeat(n, fx){
    if (n == 0){
        return [];
    }else{
        return [fx()].concat(repeat(n-1, fx));
    }
}

var tf_eq = make_factor(
        function(x, y){
            console.log([x, y]);
            if (x == y){
                console.log("eq!");
                return Math.log(0.05);
            }else{
                console.log("NOT eq!");
                return Math.log(1.0);
            }
        }
        );

console.log("-----------------");
function my_distribution(){
    var sample = repeat(2, flip);
    tf_eq(sample[0], sample[1]);
    return sample; 
}
var results = mh_query(100, 1, my_distribution);
//console.log("============================");
//console.log(stores[0]);
//console.log("============================");
//console.log(stores[1]);
console.log(results);


