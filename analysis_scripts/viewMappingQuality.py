import sys, subprocess
from collections import defaultdict

bashCommand = "samtools view " + sys.argv[1]
try:
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
except:
        print "load samtools module"
        sys.exit()

output, error = process.communicate()

lines = output.strip().split('\n')
myDict = defaultdict(int)
for line in lines:
        ll = line.split('\t')
	myDict[int(ll[4])] += 1


with open("mapq/" + sys.argv[1].replace('.bam','').replace('bamfiles/', '') + '.mapq', 'w') as outfile:
	outfile.write("mapq\tcount\n")
	for key in sorted(myDict.keys()):
		outfile.write(str(key) + '\t' + str(myDict[key]) + '\n')

