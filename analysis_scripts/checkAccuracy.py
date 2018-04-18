import sys
import subprocess
from collections import defaultdict



infile = open(sys.argv[1], 'r')
in2file = open(sys.argv[2], 'r')

# NC_000913.3-618880    83    NC_000913.3    892918    99    150=    =    892867    -201    CTGTGCCGGTCTTTTTATCGAAAGAGGTTGTACAAAATTATGACATCGCTGGTCGTTCCTGGTCTGGATACGCTGCGTCAATGGCTCGATGACCTGGGGATGAGTTTTTTTGAATGTGATAACTGTCAGGCTCTGCATCTGCCCCATATG    GGCGGGGGGGCGGGGCC8GCCCGGCGGGGCGCGGGCG8CGCGCCGGGCGJGGGCCGGCGCGGGG=G=CGGG=GGGGGGGCGGGGGJGGJJJJJJGJ1GGJCJGGJJJJG8JGGJJJJJJCGJJJJGJJGJGJJJGJJGG1GGGGGGCCCC
# NC_000913.3-618880    163    NC_000913.3    892867    99    14=1X135=    =    892918    201    TAGTCGCAATCACAGTACTGATCATGGTTTTGCCTGCGCTTTTTGCGTAAGCTGTGCCGGTCTTTTTATCGAAAGAGGTTGTACAAAATTATGACATCGCTGGTCGTTCCTGGTCTGGATACGCTGCGTCAATGGCTCGATGACCTGGGG    CCCGGGGGGGGGGJ(JJJGGGJJJJGGGJJGJGJJGJ8JGJJ1JJ8JGJJGCJGJGGCJCJCJCJGGGGC=C=CCGCG8=GGCGGCCGCCGCCGGCCG=8CCCCCJ=GGCGCGGGG=81GCGGGCCGGCGGC=GGGCGGGGG=GGGCG1=

#  ##ART_Illumina    read_length    150
#  @CM    ./art_illumina -ss HS25 -sam -i /fslhome/penrodce/compute/model/miseq_ecoli/ecoli.fasta -p -l 150 -f 200 -m 200 -s 10 -o /fslhome/penrodce/compute/model/ecoli_sim_200x -rs 1521495875
#  @SQ    NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome    4641652
#  ##Header End
#  >NC_000913.3    NC_000913.3-6188800/1-1    3731838    -
#  GCCGCAGGCGCTGGGTGCGCAGGCGACAGAGCCAGAACGTCAGGTGGTCGCCATGTGCGGCGATGGCGGTTTTAGCATGTTGATGGGCGATTTCCTCTCAGTAGTGCAGATGAAACTGCCAGTGAAAATTGTCGTCTTTAACAACAGCGT
#  GCCGCAGGCGCTGGGTGCGCAGGCGACAGAGCCAGAACGTCAGGTGGTCGCCATGTGCGGCGATGGCGGTTTTAGCATGTTTATGGGCGATTTCCTCTCAGTAGTGCAGATGAAACTGCCAGTGAAAATTGTCGTCTTTAACAACAGCGT
#  >NC_000913.3    NC_000913.3-6188798/1-1    2411070    +
#  CATTTTGATTTCACTATTGCCCGCAACTTTGACCAGCCATGTTCCCATTACGGTGGCAACACCGGTACGGACCAAACCATCGCCAATAATAAACAAGGCGGCAATCAGGACAACGTTAGGATCAGAAAAGCCGGAAAATACTTCTGGGAC
#  CATTTTGATTTCACTATTGCCCGCAACTTTGACCAGCCATGTTCCCATTACGGTGGCAACACCGGTACGGACCAAACCATCGCCAATAATAAACAAGGCGGCAATCAGGACAACGTTAGGATCAGAAAAGCCGGAAAATACTTCTGGGAC

readDict = defaultdict(list)
for line in infile:
	if line[0] == '>':
		ll = line.strip().split('\t')
		name = ll[1].split('/')[0]
		readDict[name].append(int(ll[2]))

for line in in2file:
	if line[0] == '>':
		ll = line.strip().split('\t')
		name = ll[1].split('/')[0]
		readDict[name].append(int(ll[2]))



infile.close()

outfile = open('accuracy.tsv', 'w')
outfile.write("file\taccuracy\n")
for file in sys.argv[3:]:
	bashCommand = "samtools view " + file
	process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
	
	output, error = process.communicate()
	
	lines = output.strip().split('\n')
	
	correct = 0
	total = 0
	for line in lines:
		if line[0] != '#':
			ll = line.strip().split('\t')
			if int(ll[4]) == 42:
				total += 1
				if ll[0] not in readDict:
					print "NONE"
				valList = readDict[ll[0]]
				value = int(ll[3])
				for val in valList:
					if val - 50 < value < val + 50:
						correct += 1
						break
		
	
	outfile.write(file + '\t' + str(correct/float(total)) + '\n')

outfile.close()
