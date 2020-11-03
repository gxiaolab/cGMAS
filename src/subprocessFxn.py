# 2/7/2017
# from Ashley Cass

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
