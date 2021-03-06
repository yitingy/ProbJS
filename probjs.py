import sys
import os

filenames = os.listdir('random')
filenames = filter(lambda f: f.endswith(".js"), filenames)
filenames = map(lambda f: 'random/'+f, filenames)
filenames += ['math-functions.js', 'stacktrace.js', 'xrp.js', "xrp-preamble.js", 'factor.js', 'mcmc-preamble.js'];
infile = sys.argv[1]
outfile = "temp.js"

fh = open(outfile, 'w')
for filename in filenames:
    fh.write(open(filename).read())

fh.write(open(infile).read())

fh.close()


