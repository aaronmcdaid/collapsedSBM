.PHONY: gitstatus.txt help clean
BITS=
#BITS=-m32
#BITS=-m64

all: acp

clean:
	-rm tags acp *.o

tags:
	ctags *.[ch]pp


					#-Wclobbered   -Wempty-body   \ -Wignored-qualifiers  -Woverride-init   \ -Wtype-limits   -Wunused-but-set-parameter 
# I'm including most of the -Wextra flags, but I want rid of the enum-in-conditional warning from boost
CFLAGS=       \
          -Wmissing-field-initializers   \
          -Wsign-compare   \
          -Wuninitialized   \
          -Wunused-parameter    \
          -Wunused             \
          -Wall -Wformat -Werror -I./boost_1_41_0

#CXXFLAGS= ${BITS}     -g
LDFLAGS+= -lstdc++ -lrt
CXXFLAGS= ${BITS} -O3     -g ${CFLAGS} # -DNDEBUG
#CXXFLAGS=              -O2                 

acp: shmGraphRaw.o Range.o
