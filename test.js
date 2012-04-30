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


