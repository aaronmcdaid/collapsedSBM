#!/usr/bin/python
import collections
import pylab
import random
import networkx as nx

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
	print "N=", N
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
	return sorted(all_data, key=lambda a: a[1:3])

def main():
	all_data = process_z_unswitched()
	for a in all_data:
		print a

def plot(all_data):
	#G=nx.read_weighted_edgelist('../../edge_list.txt', nodetype=int, delimiter=' ', create_using=nx.DiGraph())
	#G=nx.read_weighted_edgelist('../../edge_list.txt', nodetype=int, delimiter='\t', create_using=nx.DiGraph())
	G=nx.read_weighted_edgelist('../../edge_list.txt', nodetype=int, delimiter=' ', create_using=nx.DiGraph())
	N=G.order()
	node_set = sorted(G.nodes())
	print node_set
	assert N == len(node_set)
	print "N=",N
	assert node_set == sorted([node_set[a[0]] for a in all_data])
	mat=nx.convert.to_numpy_matrix(G, nodelist=[node_set[a[0]] for a in all_data])
	#mat=[[random.randint(0,1) for i in range(N)] for j in range(N)]
	pylab.cla(); pylab.imshow(pylab.asarray(1-mat), cmap=pylab.cm.gray, interpolation='nearest')

	current_cluster = all_data[0][1]
	for i, x in enumerate(all_data):
		if i==0:
			continue
		if current_cluster != x[1]:
			print x[1]
			current_cluster = x[1]
			pylab.axhline(y=i-0.55,color='red',linewidth=2)
			pylab.axvline(x=i-0.45,color='red',linewidth=2)


if __name__ == "__main__":
    main()
