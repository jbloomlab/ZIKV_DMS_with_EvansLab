import time
from Bio import SeqIO
import os

def main():
	start_time = time.time()
	localtime = time.asctime()
	print "Start time: %s" % localtime

	outfile_name = 'E_var_muts.txt'
	if os.path.isfile(outfile_name):
		os.remove(outfile_name)
	outfile = open(outfile_name, 'w')

	refseq_file = 'E.fa'
	refseq = SeqIO.parse(open(refseq_file), 'fasta')
	for record in refseq:
		refseq = str(record.seq)

	assert len(refseq) % 3 == 0, "refseq does not specify a sequence of codons since its length of %d nucleotides is not a multiple of 3" % (len(refseq))
	sites = [i + 1 for i in range(len(refseq) // 3)]
	poss = [i + 1 for i in range(len(refseq))]
	refseq_chars = {}
	refseq_pos = {}
	for site in sites:
		refseq_chars[site] = refseq[3 * site - 3 : 3 * site]


	muts_file = 'E_clones.fa'
	muts_seqs = SeqIO.parse(open(muts_file), 'fasta')
	for mut in muts_seqs:
		mutname = str(mut.name)
		mutseq = str(mut.seq)
		muts = []
		mut_chars = {}
		
		for site in sites:
			mut_chars[site] = mutseq[3 * site - 3 : 3 * site]
			mutation = mut_chars[site]
			if mutation != refseq_chars[site]:
				if '-' in mutation:
					for i in range(3):
						if mutation[i] == '-':
							nt_pos = site*3 - (2 - i) 
							muts.append('del%s%d' % (refseq_chars[site][i], nt_pos))
								
				else:
					muts.append(refseq_chars[site]+str('%0.3d' % site)+ mutation) 
		if len(muts) == 0:
			muts.append('no_mutations')
		outfile.write(mutname + ': '+','.join(str(m) for m in muts)+'\n')

	outfile.close()
















	interval = float((time.time() - start_time)/60)
	time_taken = float(interval)
	unit = ''
	if interval > 60:
		if interval > 1440:
			time_taken = interval/1440
			unit = 'days'
		else:
			time_taken = interval/60
			unit = 'hours'
	else: 
		tame_taken = interval
		unit = 'minutes'

	print "\n---------------------------------------- \nthe program took %.2f %s to run \n----------------------------------------\n" % (time_taken, unit)





if __name__ == '__main__':
	main() #run the script