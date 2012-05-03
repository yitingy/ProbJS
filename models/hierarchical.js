var colors = ["blue", "green", "red"];

var samples = mh_query(
  200, 100,
  function(){

    var phi = dirichlet([1, 1, 1]);
    var alpha = .1;
    var prototype = [alpha * w for (w in phi)];

    var bag_prototype = mem( function(obs_name){ dirichlet(prototype); } );

    var obs_bag =
      mem( function(obs_name){
             uniform_draw([bag1, bag2, bag3]); } );

    var draw_marble =
      mem( function(obs_name){
             multinomial(colors, bag_prototype(obs_bag(obs_name)));} );

    condition(
      [draw_marble( "obs1") == "red",
       draw_marble( "obs2") == "red",
       draw_marble( "obs3") == "blue",
       draw_marble( "obs4") == "blue",
       draw_marble( "obs5") == "red",
       draw_marble( "obs6") == "blue"]);

    return [ (obs_bag("obs1") == obs_bag("obs2")),
             (obs_bag("obs1") == obs_bag("obs3")) ];
    
  });