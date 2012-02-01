.PHONY: gitstatus.txt help clean
BITS=
#BITS=-m32
#BITS=-m64

MAIN=sbm

all: tags ${MAIN}

clean:
	-rm tags ${MAIN} *.o */*.o

tags:
	ctags *.[ch]pp


					#-Wclobbered   -Wempty-body   \ -Wignored-qualifiers  -Woverride-init   \ -Wtype-limits   -Wunused-but-set-parameter 
CXXFLAGS=       \
          -Wmissing-field-initializers   \
          -Wsign-compare   \
          -Wuninitialized   \
          -Wunused-parameter    \
          -Wunused             \
          -Wnon-virtual-dtor            \
          -Wall -Wformat -Werror -Wextra #-Wconversion # scf.cpp doesn't like -Wconversion

CXX=/home/aaronmcdaid/bin/g++-4.6
CC=/home/aaronmcdaid/bin/g++-4.6
#CXXFLAGS= ${BITS}     -g
LDFLAGS+= -lrt
LDFLAGS+= -lgsl -lgslcblas
CXXFLAGS:= ${BITS} -O3        ${CXXFLAGS} -std=gnu++0x # -DNDEBUG
#CXXFLAGS+= -p -pg
#CXXFLAGS=              -O2                 

#${MAIN}: CXXFLAGS += -DNDEBUG
${MAIN}: gitstatus.o Range.o sbm_state.o cmdline.o scf.o graph/weights.o graph/strings.o graph/loading.o graph/network.o graph/saving.o graph/bloom.o format_flag_stack/format_flag_stack.o

#lineGraph: lineGraph.o shmGraphRaw.o Range.o
gitstatus.txt: 
	{ git log | head -n 1 ; git status ; } | head -n 20 | sed -re 's/"/\\"/g ; s/^/"/g; s/$$/\\n"/g; ' > gitstatus.txt
gitstatus.o: comment.txt  gitstatus.txt

cmdline.c.FORCE:
	# remake cmdline.c . But it's OK unless you change the .ggo file. You'll need gengetopt(1) to be able to run this.
	gengetopt  --unamed-opts < cmdline.ggo
