// From http://baagoe.com/en/RandomMusings/javascript/
function Alea() {
  return (function(args) {
    // Johannes Baagøe <baagoe@baagoe.com>, 2010
    var s0 = 0;
    var s1 = 0;
    var s2 = 0;
    var c = 1;

    if (args.length == 0) {
      args = [+new Date];
    }
    var mash = Mash();
    s0 = mash(' ');
    s1 = mash(' ');
    s2 = mash(' ');

    for (var i = 0; i < args.length; i++) {
      s0 -= mash(args[i]);
      if (s0 < 0) {
        s0 += 1;
      }
      s1 -= mash(args[i]);
      if (s1 < 0) {
        s1 += 1;
      }
      s2 -= mash(args[i]);
      if (s2 < 0) {
        s2 += 1;
      }
    }
    mash = null;

    var random = function() {
      var t = 2091639 * s0 + c * 2.3283064365386963e-10; // 2^-32
      s0 = s1;
      s1 = s2;
      return s2 = t - (c = t | 0);
    };
    random.uint32 = function() {
      return random() * 0x100000000; // 2^32
    };
    random.fract53 = function() {
      return random() + 
        (random() * 0x200000 | 0) * 1.1102230246251565e-16; // 2^-53
    };
    random.version = 'Alea 0.9';
    random.args = args;
    return random;

  } (Array.prototype.slice.call(arguments)));
};
// From http://baagoe.com/en/RandomMusings/javascript/
function KISS07() {
  return (function(args) {
    // George Marsaglia, 2007-06-23
    //http://groups.google.com/group/comp.lang.fortran/msg/6edb8ad6ec5421a5
    var x = 123456789;
    var y = 362436069;
    var z =  21288629;
    var w =  14921776;
    var c = 0;

    if (args.length == 0) {
      args = [+new Date];
    }
    var mash = Mash();
    for (var i = 0; i < args.length; i++) {
      x ^= mash(args[i]) * 0x100000000; // 2^32
      y ^= mash(args[i]) * 0x100000000;
      z ^= mash(args[i]) * 0x100000000;
      w ^= mash(args[i]) * 0x100000000;
    }
    if (y === 0) {
      y = 1;
    }
    c ^= z >>> 31;
    z &= 0x7fffffff;
    if ((z % 7559) === 0) {
      z++;
    }
    w &= 0x7fffffff;
    if ((w % 7559) === 0) {
      w++;
    }
    mash = null;

    var uint32 = function() {
      var t;

      x += 545925293;
      x >>>= 0;

      y ^= y << 13;
      y ^= y >>> 17;
      y ^= y << 5;

      t = z + w + c;
      z = w;
      c = t >>> 31;
      w = t & 0x7fffffff;

      return x + y + w >>> 0;
    };

    var random = function() {
      return uint32() * 2.3283064365386963e-10; // 2^-32
    };
    random.uint32 = uint32;
    random.fract53 = function() {
      return random() +
        (uint32() & 0x1fffff) * 1.1102230246251565e-16; // 2^-53
    };
    random.args = args;
    random.version = 'KISS07 0.9';

    return random;
  } (Array.prototype.slice.call(arguments)));
};
// From http://baagoe.com/en/RandomMusings/javascript/
function Kybos() {
  return (function(args) {
    // Johannes Baagøe <baagoe@baagoe.com>, 2010
    var s0 = 0;
    var s1 = 0;
    var s2 = 0;
    var c = 1;
    var s = [];
    var k = 0;

    var mash = Mash();
    var s0 = mash(' ');
    var s1 = mash(' ');
    var s2 = mash(' ');
    for (var j = 0; j < 8; j++) {
      s[j] = mash(' ');
    }

    if (args.length == 0) {
      args = [+new Date];
    }
    for (var i = 0; i < args.length; i++) {
      s0 -= mash(args[i]);
      if (s0 < 0) {
        s0 += 1;
      }
      s1 -= mash(args[i]);
      if (s1 < 0) {
        s1 += 1;
      }
      s2 -= mash(args[i]);
      if (s2 < 0) {
        s2 += 1;
      }
      for (var j = 0; j < 8; j++) {
        s[j] -= mash(args[i]);
        if (s[j] < 0) {
          s[j] += 1;
        }
      }
    }

    var random = function() {
      var a = 2091639;
      k = s[k] * 8 | 0;
      var r = s[k];
      var t = a * s0 + c * 2.3283064365386963e-10; // 2^-32
      s0 = s1;
      s1 = s2;
      s2 = t - (c = t | 0);
      s[k] -= s2;
      if (s[k] < 0) {
        s[k] += 1;
      }
      return r;
    };
    random.uint32 = function() {
      return random() * 0x100000000; // 2^32
    };
    random.fract53 = function() {
      return random() +
        (random() * 0x200000 | 0) * 1.1102230246251565e-16; // 2^-53
    };
    random.addNoise = function() {
      for (var i = arguments.length - 1; i >= 0; i--) {
        for (j = 0; j < 8; j++) {
          s[j] -= mash(arguments[i]);
          if (s[j] < 0) {
            s[j] += 1;
          }
        }
      }
    };
    random.version = 'Kybos 0.9';
    random.args = args;
    return random;

  } (Array.prototype.slice.call(arguments)));
};
// From http://baagoe.com/en/RandomMusings/javascript/
function LFib() {
  return (function(args) {
    // Johannes Baagøe <baagoe@baagoe.com>, 2010
    var k0 = 255,
        k1 = 52,
        k2 = 0;
    var s = [];

    var mash = Mash();
    if (args.length === 0) {
      args = [+new Date()];
    }
    for (var j = 0; j < 256; j++) {
      s[j] = mash(' ');
      s[j] -= mash(' ') * 4.76837158203125e-7; // 2^-21
      if (s[j] < 0) {
        s[j] += 1;
      }
    }
    for (var i = 0; i < args.length; i++) {
      for (var j = 0; j < 256; j++) {
        s[j] -= mash(args[i]);
        s[j] -= mash(args[i]) * 4.76837158203125e-7; // 2^-21
        if (s[j] < 0) {
          s[j] += 1;
        }
      }
    }
    mash = null;

    var random = function() {
      k0 = (k0 + 1) & 255;
      k1 = (k1 + 1) & 255;
      k2 = (k2 + 1) & 255;

      var x = s[k0] - s[k1];
      if (x < 0.0) {
        x += 1.0;
      }
      x -= s[k2];
      if (x < 0.0) {
        x += 1.0;
      }
      return s[k0] = x;
    }

    random.uint32 = function() {
      return random() * 0x100000000 >>> 0; // 2^32
    };
    random.fract53 = random;
    random.version = 'LFib 0.9';
    random.args = args;

    return random;
  } (Array.prototype.slice.call(arguments)));
};
// From http://baagoe.com/en/RandomMusings/javascript/
function LFIB4() {
  return(function(args) {
    // George Marsaglia's LFIB4,
    //http://groups.google.com/group/sci.crypt/msg/eb4ddde782b17051
    var k0 = 0,
        k1 = 58,
        k2 = 119,
        k3 = 178;

    var s = [];

    var mash = Mash();
    if (args.length === 0) {
      args = [+new Date()];
    }
    for (var j = 0; j < 256; j++) {
      s[j] = mash(' ');
      s[j] -= mash(' ') * 4.76837158203125e-7; // 2^-21
      if (s[j] < 0) {
        s[j] += 1;
      }
    }
    for (var i = 0; i < args.length; i++) {
      for (var j = 0; j < 256; j++) {
        s[j] -= mash(args[i]);
        s[j] -= mash(args[i]) * 4.76837158203125e-7; // 2^-21
        if (s[j] < 0) {
          s[j] += 1;
        }
      }
    }
    mash = null;

    var random = function() {
      var x;

      k0 = (k0 + 1) & 255;
      k1 = (k1 + 1) & 255;
      k2 = (k2 + 1) & 255;
      k3 = (k3 + 1) & 255;

      x = s[k0] - s[k1];
      if (x < 0) {
        x += 1;
      }
      x -= s[k2];
      if (x < 0) {
        x += 1;
      }
      x -= s[k3];
      if (x < 0) {
        x += 1;
      }

      return s[k0] = x;
    }

    random.uint32 = function() {
      return random() * 0x100000000 >>> 0; // 2^32
    };
    random.fract53 = random;
    random.version = 'LFIB4 0.9';
    random.args = args;

    return random;
  } (Array.prototype.slice.call(arguments)));
};
// From http://baagoe.com/en/RandomMusings/javascript/
// Johannes Baagøe <baagoe@baagoe.com>, 2010
function Mash() {
  var n = 0xefc8249d;

  var mash = function(data) {
    data = data.toString();
    for (var i = 0; i < data.length; i++) {
      n += data.charCodeAt(i);
      var h = 0.02519603282416938 * n;
      n = h >>> 0;
      h -= n;
      h *= n;
      n = h >>> 0;
      h -= n;
      n += h * 0x100000000; // 2^32
    }
    return (n >>> 0) * 2.3283064365386963e-10; // 2^-32
  };

  mash.version = 'Mash 0.9';
  return mash;
}

// From http://baagoe.com/en/RandomMusings/javascript/
function MRG32k3a() {
  return (function(args) {
    // Copyright (c) 1998, 2002 Pierre L'Ecuyer, DIRO, Université de Montréal.
    // http://www.iro.umontreal.ca/~lecuyer/
    var m1 = 4294967087;
    var m2 = 4294944443;
    var s10 = 12345,
        s11 = 12345,
        s12 = 123,
        s20 = 12345,
        s21 = 12345,
        s22 = 123;

    if (args.length === 0) {
      args = [+new Date()];
    }
    var mash = Mash();
    for (var i = 0; i < args.length; i++) {
      s10 += mash(args[i]) * 0x100000000; // 2 ^ 32
      s11 += mash(args[i]) * 0x100000000;
      s12 += mash(args[i]) * 0x100000000;
      s20 += mash(args[i]) * 0x100000000;
      s21 += mash(args[i]) * 0x100000000;
      s22 += mash(args[i]) * 0x100000000;
    }
    s10 %= m1;
    s11 %= m1;
    s12 %= m1;
    s20 %= m2;
    s21 %= m2;
    s22 %= m2;
    mash = null;

    var uint32 = function() {
      var m1 = 4294967087;
      var m2 = 4294944443;
      var a12 = 1403580;
      var a13n = 810728;
      var a21 = 527612;
      var a23n = 1370589;

      var k, p1, p2;

      /* Component 1 */
      p1 = a12 * s11 - a13n * s10;
      k = p1 / m1 | 0;
      p1 -= k * m1;
      if (p1 < 0) p1 += m1;
      s10 = s11;
      s11 = s12;
      s12 = p1;

      /* Component 2 */
      p2 = a21 * s22 - a23n * s20;
      k = p2 / m2 | 0;
      p2 -= k * m2;
      if (p2 < 0) p2 += m2;
      s20 = s21;
      s21 = s22;
      s22 = p2;

      /* Combination */
      if (p1 <= p2) return p1 - p2 + m1;
      else return p1 - p2;
    };

    var random = function() {
      return uint32() * 2.3283064365386963e-10; // 2^-32
    };
    random.uint32 = uint32;
    random.fract53 = function() {
      return random() +
        (uint32() & 0x1fffff) * 1.1102230246251565e-16; // 2^-53
    };
    random.version = 'MRG32k3a 0.9';
    random.args = args;

    return random;
  } (Array.prototype.slice.call(arguments)));
};
// From http://baagoe.com/en/RandomMusings/javascript/
function Xorshift03() {
  return (function(args) {
    // George Marsaglia, 13 May 2003
    // http://groups.google.com/group/comp.lang.c/msg/e3c4ea1169e463ae
    var x = 123456789,
        y = 362436069,
        z = 521288629,
        w = 88675123,
        v = 886756453;

    if (args.length == 0) {
      args = [+new Date];
    }
    var mash = Mash();
    for (var i = 0; i < args.length; i++) {
      x ^= mash(args[i]) * 0x100000000; // 2^32
      y ^= mash(args[i]) * 0x100000000;
      z ^= mash(args[i]) * 0x100000000;
      v ^= mash(args[i]) * 0x100000000;
      w ^= mash(args[i]) * 0x100000000;
    }
    mash = null;

    var uint32 = function() {
      var t = (x ^ (x >>> 7)) >>> 0;
      x = y;
      y = z;
      z = w;
      w = v;
      v = (v ^ (v << 6)) ^ (t ^ (t << 13)) >>> 0;
      return ((y + y + 1) * v) >>> 0;
    }

    var random = function() {
      return uint32() * 2.3283064365386963e-10; // 2^-32
    };
    random.uint32 = uint32;
    random.fract53 = function() {
      return random() +
        (uint32() & 0x1fffff) * 1.1102230246251565e-16; // 2^-53
    };
    random.version = 'Xorshift03 0.9';
    random.args = args;
    return random;

  } (Array.prototype.slice.call(arguments)));
};
// Set the random number generator
var random = new MRG32k3a();
var intRandom = random.uint32;

// Use these functions from the scheme2js runtime.js if available
var sc_list2vector = (typeof sc_list2vector == "function") ? sc_list2vector : function(x) { return x };
var sc_vector2list = (typeof sc_vector2list == "function") ? sc_vector2list : function(x) { return x };

function random_integer(n)
{ 
    return intRandom() % n; 
}

function random_real()
{ 
    return random(); 
}

function seed_rng(seed)
{
    random = new MRG32k3a(seed);
    intRandom = random.uint32;
}

// Draw sample from Poisson distribution
// Knuth TAOCP 2 (roughly optimal)
function sample_poisson(mu)
{
    var k = 0;

    while(mu > 10)
    {
        var m = 7/8*mu;
        var x = Math.sample_gamma(m);

        if(x > mu) return k + sample_binomial(mu/x, m-1);
        else{ mu -= x; k += m; }
    }

    var emu = Math.exp(-mu);
    var p = 1;
    do{ p *= random(); k++; } while(p > emu);

    return k-1;
}

// Poisson probability distribution function via iterative expansion
function poisson_pdf(k, mu)
{
    return Math.exp(k * Math.log(mu) - mu - lnfact(k));
}

// Draw sample from a Gamma distribution
// Marsagli and Tsang '00 (roughly optimal)
function sample_gamma(a,b)
{
    if(a < 1) return sample_gamma(1+a,b) * Math.pow(random(), 1/a);

    var x,v,u;
    var d = a-1/3;
    var c = 1/Math.sqrt(9*d);

    while(true)
    {
        do{x = sample_gaussian(0,1);  v = 1+c*x;} while(v <= 0);

        v=v*v*v;
        u=random();

        if((u < 1 - .331*x*x*x*x) || (Math.log(u) < .5*x*x + d*(1 - v + Math.log(v)))) return b*d*v;
    }
}

// Evaluate gamma pdf
function gamma_pdf(x,a,b)
{
    if(x<0) return 0;
    if(x==0) return a==1 ? 1/b : 0;
    if(a==1) return Math.exp(-x/b)/b;
    
    return Math.exp((a - 1)*Math.log(x/b) - x/b - log_gamma(a))/b;
}

// Evaluate log gammma pdf
function gamma_lnpdf(x,a,b)
{
    return (1 - a)*Math.log(x) - x/b - log_gamma(a) - a*Math.log(b);
}

// Draw a sample from a Binomial distribution
// Knuth TAOCP 2 (could be improved, i.e. via Kachitvichyanukul & Schmeiser)
function sample_binomial(p,n)
{
    var k = 0;
    var N = 10;

    var a, b;
    while(n > N)
    {
        a = 1 + n/2;
        b = 1 + n-a;

        var x = sample_beta(a,b);

        if(x >= p){ n = a-1; p /= x; }
        else{ k += a; n = b - 1; p = (p-x) / (1-x); }
    }

    var u;
    for(i=0; i<n; i++)
    {
        u = random();
        if(u<p) k++;
    }

    return k;
}

// Binomial probability distribution function via Normal approximation
// Peizer & Pratt 1968, JASA 63: 1416-1456 (may not be optimal...)
function binomial_pdf(k, p, n)
{
    var inv2 = 1/2, inv3 = 1/3, inv6 = 1/6;

    if (k >= n) return 1;

    var q = 1 - p;
    var s = k + inv2;
    var t = n - k - inv2;
    var d1 = s + inv6 - (n + inv3) * p;
    var d2 = q /(s+inv2) - p/(t+inv2) + (q-inv2)/(n+1);

    d2 = d1 + 0.02 * d2;

    var num = 1 + q * g(s/(n*p)) + p * g(t/(n*q));
    var den = (n + inv6) * p * q;
    var z = num / den;

    z = d2 * Math.sqrt(z);
    z = normal_cdf(z);

    return z;
}

// Draw a sample from a Beta distribution
// Knuth TAOCP 2 (roughly optimal)
function sample_beta(a, b)
{
    var x = sample_gamma(a, 1);
    return x / (x + sample_gamma(b, 1));
}

// Draw a sample from a Gaussian distribution
// Leva '92 (could be improved, i.e. via Ziggurat method)
function sample_gaussian(mu,sigma)
{
    var u, v, x, y, q;

    do
    {
        u = 1 - random();
        v = 1.7156 * (random() - .5);
        x = u - 0.449871;
        y = Math.abs(v) + 0.386595;
        q = x*x + y*(0.196*y - 0.25472*x);
    }
    while(q >= 0.27597 && (q > 0.27846 || v*v > -4 * u * u * Math.log(u)))

    return mu + sigma*v/u;
}

// Evaluate the gaussian distribution
function gaussian_pdf(x,mu,sigma)
{
    x-=mu;
    var asigma = Math.abs(sigma);
    var u = x/asigma;
    return (1/ Math.sqrt(2*Math.PI) * asigma) * Math.exp(-u*u/2);  
}

// Evaluate the log gaussian distribution
function gaussian_lnpdf(x,mu,sigma)
{
    return -.5*(1.8378770664093453 + Math.log(sigma) + (x - mu)*(x - mu)/sigma);
}

// Draw a sample from a Dirichlet distribution
// Law & Kelton (roughly optimal)
// TODO: may need to match function signature for Ikarus compatibility
// TODO: handle underflow in normalization
function sample_dirichlet(alpha)
{
    alpha = sc_list2vector(alpha);
    var theta = new Array(alpha.length);
    var sum = 0;

    for(i=0; i<alpha.length; i++){ theta[i] = sample_gamma(alpha[i],1); sum += theta[i]; }
    for(i=0; i<alpha.length; i++) theta[i] /= sum;
    
    return sc_vector2list(theta);
}

// Evaluate the logarithm of the Dirichlet distribution
function dirichlet_lnpdf(theta, alpha)
{
    alpha = sc_list2vector(alpha);
    theta = sc_list2vector(theta);
    var logp = log_gamma(sum(alpha));
    
    for(i=0; i<alpha.length; i++) logp += (alpha[i] - 1)*Math.log(theta[i]);
    for(i=0; i<alpha.length; i++) logp -= log_gamma(alpha[i]);

    return logp;      
}

// Draw a sample from a Student's t-distribution
// Marsaglia '80
function sample_tdist(nu)
{
    if(nu <= 2) return sample_gaussian(0,1) / sqrt( 2 * sample_gamma(nu/2, 1) / nu);

    var a,b,c,t;
    do
    {
        a = sample_gaussian(0,1);
        b = -1 / (nu/2 - 1) * log1p(-random());
        c = a*a/(nu - 2);
    }
    while(1-c < 0 || Math.exp(-b-c) > (1-c));

    return a / Math.sqrt((1-c/nu) * (1-c));
}

// Evaluate t-distribution
function tdist_pdf(x,nu)
{
    var a = log_gamma(nu/2);
    var b = log_gamma((nu+1)/2);
    
    return Math.exp(b-a)/Math.sqrt(Math.PI*nu) * Math.pow(1 + x*x/nu, -(nu+1)/2);
}

// Draw a sample from a generalized t-distribution
function sample_generalized_tdist(nu,mu,sigma_squared)
{
    return sample_tdist(nu)*Math.sqrt(sigma_squared) + mu;
}

// Return the log of a sum of exponentials, to minimize under/overflow
function logsumexp(v)
{
    v = sc_list2vector(v);
    var t=0,
        val;

    for(i=0;i<v.length;i++)
    {
        var abs=Math.abs(v[i]);        
        if(abs>t){ t=abs; val=v[i]; }                          
    }

    var sum=0;
    for(i=0;i<v.length;i++) {
      sum += Math.exp(v[i]-val);
    }

    return Math.log(sum) + val;
}

// Evaluate the log of gamma(x)
// Lancsoz approximation from Numerical Recipes in C
function log_gamma(xx)
{
    var cof = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5]; 

    var x = xx - 1.0;
    var tmp = x + 5.5; tmp -= (x + 0.5)*Math.log(tmp);
    var ser=1.000000000190015;
    for (j=0;j<=5;j++){ x++; ser += cof[j]/x; }
    return -tmp+Math.log(2.5066282746310005*ser);
}

// Calculate the sum of elements in a vector
// N.B.: this doesn't get used in compiled Church->JS code
// so we don't need to use sc_list2vector / sc_vector2list
function sum(v)
{
    var sum=0;
    for(i=0;i<v.length;i++) sum += v[i];
    return sum;
}

// Calculate the mean of elements in a vector
// N.B.: this doesn't get used in compiled Church->JS code
// so we don't need to use sc_list2vector / sc_vector2list
function mean(v)
{
    return sum(v)/v.length;
}

// Normalize a vector
function normalize(v)
{
    v = sc_list2vector(v);
    var s=0;
    for(i=0;i<v.length;i++) s += v[i]*v[i];
    s = Math.sqrt(s);
    for(i=0;i<v.length;i++) v[i] /= s;
    return sc_vector2list(v);
}

// Returns log(1 + x) in a numerically stable way
function log1p(x)
{
    var ret = 0;
    var n = 50; // degree of precision

    if(x <= -1) return Number.NEGATIVE_INFINITY;
    if(x < 0 || x > 1) return Math.log(1+x);

    for(i=1; i<n; i++)
        if ((i % 2) === 0) ret -= Math.pow(x,i)/i;
        else ret += Math.pow(x,i)/i;

    return ret;
}

// factorial(x)
function fact(x)
{
    var t=1;
    while(x>1) t*=x--;
    return t;
}

// ln(x!) by Stirling's formula
// [Knuth I: p111]
function lnfact(x)
{
    if (x < 1) x = 1;

    if (x < 12) return Math.log(fact(Math.round(x)));

    var invx = 1 / x;
    var invx2 = invx * invx;
    var invx3 = invx2 * invx;
    var invx5 = invx3 * invx2;
    var invx7 = invx5 * invx2;

    var sum = ((x + 0.5) * Math.log(x)) - x;
    sum += Math.log(2*Math.PI) / 2;
    sum += (invx / 12) - (invx3 / 360);
    sum += (invx5 / 1260) - (invx7 / 1680);

    return sum;
}

// logistic(x)
function logistic(x)
{
    return 1 / (1 + Math.exp(-x));
}

// Normal cumulative distribution function
// Abramowitz & Stegun 26.2.19
// |e(x)| < 1.5E-7
function normal_cdf(x)
{
    var d1 = 0.0498673470;
    var d2 = 0.0211410061;
    var d3 = 0.0032776263;
    var d4 = 0.0000380036;
    var d5 = 0.0000488906;
    var d6 = 0.0000053830;
    var a = Math.abs(x);
    var t;

   t = 1.0 + a*(d1+a*(d2+a*(d3+a*(d4+a*(d5+a*d6)))));

   t *= t;  t *= t;  t *= t;  t *= t;
   t = 1.0 / (t+t);

   if (x >= 0)  t = 1-t;
   return t;
}

// Peizer & Pratt 1968, JASA 63: 1416-1456
function g(x)
{
    var  switchlev = 0.1;
    var z;

    if (x == 0)  return 1;
    if (x == 1)  return 0;

    var d = 1 - x;

    if (Math.abs(d) > switchlev) return (1 - (x * x) + (2 * x * Math.log(x))) / (d * d);

    z = d / 3;
    var di = d;

    for (var i = 2; i <= 7; i++)
    {
        di *= d;
        z += (2 * di) / ((i+1) * (i+2));
    }
    return z;
}
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

var gaussian = make_stateless_xrp(
        "gaussian",
        sample_gaussian,
        function(args, val){
            return gaussian_lnpdf(val, args[0], args[1])
        },
        []
        )
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

var tf_eq_float = make_factor(
        function(x, y){
            return gaussian_lnpdf((x-y), 0, 0.01);
        }
        );

var tf_eq = make_factor(
        function(x, y){
            if (x == y){
                return Math.log(1.0);
            }else{
                return Math.log(0.05);
            }
        }
        );

console.log("-----------------");
function my_distribution(){
    var sample = repeat(2, function(){return gaussian(0, 1);});
    //var sample = repeat(2, flip);
    tf_eq_float(sample[0], sample[1]);
    return sample; 
}
var results = mh_query(100, 1, my_distribution);
//console.log("============================");
//console.log(stores[0]);
//console.log("============================");
//console.log(stores[1]);
console.log(results);


