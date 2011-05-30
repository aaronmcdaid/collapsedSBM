.PHONY: gitstatus.txt help clean
BITS=
#BITS=-m32
#BITS=-m64

MAIN=sbm

all: tags boost_1_41_0 ${MAIN}

clean:
	-rm tags ${MAIN} *.o */*.o

tags:
	ctags *.[ch]pp


					#-Wclobbered   -Wempty-body   \ -Wignored-qualifiers  -Woverride-init   \ -Wtype-limits   -Wunused-but-set-parameter 
# I'm including most of the -Wextra flags, but I want rid of the enum-in-conditional warning from boost
CXXFLAGS=       \
          -Wmissing-field-initializers   \
          -Wsign-compare   \
          -Wuninitialized   \
          -Wunused-parameter    \
          -Wunused             \
          -Wnon-virtual-dtor            \
          -Wall -Wformat -Werror -I./boost_1_41_0

boost_1_41_0:
	@echo "   " This needs Boost. It has been tested with boost 1.41 .
	@echo "   " Extract this to a folder called boost_1_41_0 . 
	@echo "   " http://sourceforge.net/projects/boost/files/boost/1.41.0/
	false

#CXXFLAGS= ${BITS}     -g
LDFLAGS+= -lstdc++ -lrt
LDFLAGS+= -lgsl -lgslcblas
CXXFLAGS:= ${BITS} -O3        ${CXXFLAGS} # -DNDEBUG
#CXXFLAGS+= -p -pg
#CXXFLAGS=              -O2                 

#${MAIN}: CXXFLAGS += -DNDEBUG
${MAIN}: gitstatus.o Range.o aaron_utils.o sbm_state.o cmdline.o scf.o graph/weights.o graph/strings.o graph/loading.o graph/network.o graph/saving.o graph/bloom.o

#lineGraph: lineGraph.o shmGraphRaw.o Range.o
gitstatus.txt: 
	{ git log | head -n 1 ; git status ; } | head -n 20 | sed -re 's/"/\\"/g ; s/^/"/g; s/$$/\\n"/g; ' > gitstatus.txt
gitstatus.o: comment.txt  gitstatus.txt

cmdline.c:
	# remake cmdline.c . But it's OK unless you change the .ggo file. You'll need gengetopt(1) to be able to run this.
	gengetopt  --unamed-opts < cmdline.ggo
