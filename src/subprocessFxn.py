# 2/7/2017
# from Ashley

# Run subprocess with input file
# cmd = ['sort', '-k1,1', '-k2,2n', infile]
# stdout, stderr = run_command(cmd)

# Run subprocess with string as input
# stdin = 'chr1\t1\t10\nchr1\t10\t20\n'
# cmd = ['sort', '-k1,1', '-k2,2n']
# stdout, stderr = run_command(cmd, stdin)

import sys
from subprocess import Popen, PIPE
def run_command(cmd, stdin=0):
	if stdin == 0:
		p = Popen(cmd, stdout=PIPE, stderr=PIPE)
		stdout, stderr = p.communicate()
	else:
		p = Popen(cmd, stdout=PIPE, stderr=PIPE, stdin=PIPE)
		stdout, stderr = p.communicate(stdin)
	if p.returncode != 0:
#		sys.stderr.write('command  failed\n')
		sys.stderr.write(stdout)
		sys.stderr.write(stderr)
		#sys.exit(1)
		stderr = 1
	return stdout, stderr

#def inf2bed():
#	#less all-er-0-dpsi-0.1-edits.txt |grep chr|grep "xon\|UTR"|cut -f1,8,9|awk -F'[\t.]' '{print $1,$4-1,$5,"exon","0",$3}' OFS='\t'|sort -u -k1,1 -k2,2n
#	stdout, stderr = subprocessFxn.run_command(['less', opts.i])
#	if stderr: print stderr; sys.exit()
#
#	stdout, stderr = subprocessFxn.run_command(['grep','chr'], stdout)
#	if stderr: print stderr; sys.exit()
#
#	# stdout, stderr = subprocessFxn.run_command(['grep','xon\|UTR'], stdout)
#	# if stderr: print stderr; sys.exit()
#
#	stdout, stderr = subprocessFxn.run_command(['cut','-f1,8,9'], stdout)
#	if stderr: print stderr; sys.exit()
#
#	stdout, stderr = subprocessFxn.run_command(['awk','-F[\t.]','{print $1,$4-1,$5,"exon","0",$3}','OFS=\t'], stdout)
#	if stderr: print stderr; sys.exit()
#
#	stdout, stderr = subprocessFxn.run_command(['sort','-u'], stdout)
#	if stderr: print stderr; sys.exit()
#
#	stdout, stderr = subprocessFxn.run_command(['sort','-k1,1','-k2,2n'], stdout)
#	if stderr: print stderr; sys.exit()
#
#	return stdout

