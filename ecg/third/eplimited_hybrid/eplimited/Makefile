# file: Makefile	G. Moody	10 April 2003
#
# make description file for compiling easytest & bxbep using gcc
#
# Install the WFDB package (http://www.physionet.org/physiotools/wfdb.shtml)
# before attempting to compile these programs.
#
# To compile these programs, just type 'make' at a command prompt while in
# this directory.  This has been tested under GNU/Linux and MS-Windows/Cygwin,
# and it should also work under MacOS/X (Darwin) and all versions of Unix.
#
# See '00README' in this directory for further information.

# CC is the name of your ANSI/ISO C compiler (might be 'cc' on some platforms).
CC = g++

qrsfilt.o:	qrsfilt-inl.h qrsdet-params.h
	$(CC) -c qrsfilt-inl.h -lm

qrsdet2.o:	qrsdet2-inl.h qrsdet-params.h
	$(CC) -c qrsdet2-inl.h
clean:
	rm -f *.o
