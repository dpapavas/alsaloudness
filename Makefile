# Set to yes, for a debug build.
DEBUG=no

# Set to yes, for an optimized build with debug output.
BENCHMARK=no

# Set to no, to disable multithreading.
THREADED=yes

# Set to yes, to enable dumping of filter response plots.
PLOT=no

# Set to yes, to build lib32 version for 64-bit systems.
LIB32=no

PREFIX = /usr/local
LIBDIR = $(PREFIX)/lib

CFLAGS = -fPIC -DPIC
CFLAGS += -Wall -Wextra -Wpedantic -Wno-unused-parameter

LDFLAGS=-lm -lasound -lfftw3f
INSTALL=/usr/bin/install

ifneq ($(DEBUG),yes)
CFLAGS += -O3
else
CFLAGS += -g
endif

ifeq ($(PLOT),yes)
CFLAGS += -DPLOT_FILTER
endif

ifeq ($(filter yes,$(DEBUG) $(BENCHMARK)),)
CFLAGS += -DNDEBUG
endif

ifeq ($(THREADED),yes)
CFLAGS += -DWITH_THREADS
LDFLAGS += -lfftw3f_threads
endif

ifeq ($(LIB32),yes)
CFLAGS += -m32
LIBDIR = $(PREFIX)/lib32
endif

all: libasound_module_pcm_loudness.so

libasound_module_pcm_loudness.so: loudness_pcm.c contours.h
	$(CC) -o $@ -shared ${CFLAGS} loudness_pcm.c ${LDFLAGS}

install: libasound_module_pcm_loudness.so
	mkdir -p $(LIBDIR)
	$(INSTALL) -m 644 libasound_module_pcm_loudness.so $(LIBDIR)/

uninstall:
	rm -f $(LIBDIR)/libasound_module_pcm_loudness.so

dist:
	if [ -e /tmp/loudness ]; then rm -rf /tmp/loudness; fi
	mkdir /tmp/loudness
	cp loudness_pcm.c contours.h Makefile /tmp/loudness
	cd /tmp; tar zcf alsaloudness.tar.gz loudness/

clean:
	rm -f libasound_module_pcm_loudness.so *~
