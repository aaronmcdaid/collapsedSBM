#!/usr/bin/Rscript
args = commandArgs(TRUE)
# print(args[1])
d = read.table(args[1], header=FALSE, sep='\t', colClasses=c('character', 'integer', 'double','double'))
N = d$V2[1]
lines = dim(d)[1]-1
stopifnot(lines %% N == 0)
timepoints = lines / N
node_names = d$V1[2:(1+N)]

all_timepoints <- vector("list", timepoints)

for(t in 0:(timepoints-1) ) {
	# print(t)
	stopifnot(all(d$V1[(t*N+2):((t+1)*N+1)] == node_names))
	# print(d$V3[(t*N+2):((t+1)*N+1)])
	# print(d$V4[(t*N+2):((t+1)*N+1)])
	# print( matrix( c(d$V3[(t*N+2):((t+1)*N+1)],d$V4[(t*N+2):((t+1)*N+1)]), nrow=N ) )
	m = matrix( c(d$V3[(t*N+2):((t+1)*N+1)],d$V4[(t*N+2):((t+1)*N+1)]), nrow=N )
	all_timepoints[[t+1]] <- m
}
length(all_timepoints)

all_data = list()
all_data.N = N
all_data.ts = all_timepoints

all_data.ts[[1]]
all_data.ts[[timepoints]]

library(MCMCpack)
procrustes(all_data.ts[[timepoints]] , all_data.ts[[1]])$X.new
