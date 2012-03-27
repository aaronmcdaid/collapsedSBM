#!/usr/bin/python
import collections

def main():
	z_unswitched = open('z.unswitched')
	header = z_unswitched.readline()
	N = int(header.strip().split(':')[1].split(',')[0])
	print N
	num_lines = 0
	counts = [collections.defaultdict(int) for n in range(N)]
	for line in z_unswitched:
		if line[0] != ' ':
			continue
		num_lines = num_lines+1
		#if num_lines > 10000:
			#break
		z = map(int,line.split('\t')[0].split())
		# print len(z), z
		assert len(z) == N
		for i, z_i in enumerate(z):
			# print i,z_i
			counts[i][z_i] = counts[i][z_i]+1
	print num_lines, counts
	print counts[67]
	print map(len, counts)


if __name__ == "__main__":
    main()
