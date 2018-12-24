.PHONY: gitstatus.txt help clean
BITS=
#BITS=-m32
#BITS=-m64

MAIN=sbm

all: ${MAIN}

clean:
	-rm ${MAIN} *.o */*.o


CXXFLAGS=       \
          -std=gnu++0x   \
          -Wmissing-field-initializers   \
          -Wsign-compare   \
          -Wuninitialized   \
          -Wunused-parameter    \
          -Wunused             \
          -Wnon-virtual-dtor            \
          -Wall -Wformat -Wextra # -Werror -Wconversion # scf.cpp doesn't like -Wconversion

CXX=g++
CC=g++ ${PROFILE}
#CXXFLAGS= ${BITS}     -g
LDFLAGS+= -lrt
LDFLAGS+= `gsl-config --libs`
CXXFLAGS:= ${BITS} -O3        ${CXXFLAGS} -std=gnu++0x # -DNDEBUG
# PROFILE= -p -pg
CXXFLAGS+= ${PROFILE}
#CXXFLAGS=              -O2                 

#${MAIN}: CXXFLAGS += -DNDEBUG
${MAIN}: gitstatus.o Range.o sbm_state.o cmdline.o scf.o graph/weights.o graph/strings.o graph/loading.o graph/network.o graph/saving.o graph/bloom.o format_flag_stack/format_flag_stack.o sbm.o
	${CXX} ${CXXFLAGS} $^ ${LDFLAGS} -o $@

#lineGraph: lineGraph.o shmGraphRaw.o Range.o
gitstatus.txt: 
	{ git log | head -n 1 ; git status ; } | head -n 20 | sed -e 's/"/\\"/g ; s/^/"/g; s/$$/\\n"/g; ' > gitstatus.txt
gitstatus.o: comment.txt  gitstatus.txt

cmdline.c.FORCE:
	# remake cmdline.c . But it's OK unless you change the .ggo file. You'll need gengetopt(1) to be able to run this.
	gengetopt  --unamed-opts < cmdline.ggo
