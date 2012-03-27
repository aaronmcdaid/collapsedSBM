#!/usr/bin/python
import collections

def average_key_in_dict(d):
	xs = [k*v for k,v in d.iteritems()]
	return (1.+sum(xs)) / sum([v for k,v in d.iteritems()])
def key_with_max_entry_in_dict(d):
	if len(d) == 1:
		return d.keys()[0]
	d2 = sorted([(v,k) for k,v in d.iteritems()])
	assert d2[-1][0] > d2[-2][0]
	return d2[-1][1]

def process_z_unswitched():
	z_unswitched = open('z.unswitched')
	header = z_unswitched.readline()
	N = int(header.strip().split(':')[1].split(',')[0])
	print N
	num_lines = 0
	counts = [collections.defaultdict(int) for n in range(N)]
	for line in z_unswitched:
		if line[0] != ' ':
			continue
		#if num_lines >= 10000:
			#break
		num_lines = num_lines+1
		z = map(int,line.split('\t')[0].split())
		# print len(z), z
		assert len(z) == N
		for i, z_i in enumerate(z):
			# print i,z_i
			counts[i][z_i] = counts[i][z_i]+1
	print num_lines
	#print counts[67]
	#print map(len, counts)
	all_data = [( n,key_with_max_entry_in_dict(counts[n]),average_key_in_dict(counts[n]),counts[n] ) for n in range(N)]
	return all_data

def main():
	all_data = process_z_unswitched()
	for a in all_data:
		print a


if __name__ == "__main__":
    main()
