CC = g++

HFILES =        CKmers.h \
        CModel.h \
        CSequences.h \
        CScore.h

CFILES =        CKmers.cpp \
        CModel.cpp \
        CSequences.cpp \
        CScore.cpp \
	#CScore.cu \
	#scoreGenome.cpp \
	#extractKmers.cpp \
        #scoreGenome.cpp \
        test.cpp

OFILES =        CKmers.o \
        CModel.o \
        CSequences.o \
        CScore.o \
	scoreReads.o \
	scoreGenome.o \
	#extractKmers.o \
        #scoreGenome.o \
        test.o

all:	${OFILES} main

main:	${OFILES}
	nvcc ${OFILES} -o main
#test:   test.o $(OFILES)
#	$(CC) -o test $(OFILES)

#test.o: CScore.h test.cpp
#	$(CC) -c test.cpp

CKmers.o:       CKmers.h CKmers.cpp
	$(CC) -c CKmers.cpp

CModel.o:       CModel.h CModel.cpp
	$(CC) -c CModel.cpp

CSequences.o:       CSequences.h CSequences.cpp
	$(CC) -c CSequences.cpp

#extractKmers.o:       extractKmers.cpp
#	$(CC) -c extractKmers.cpp

scoreGenome.o:       scoreGenome.cpp
	$(CC) -c scoreGenome.cpp

#CScore.o:	CScore.h CScore.cu
#	nvcc -c CScore.cu

CScore.o:       CScore.h CScore.cpp
	$(CC) -c CScore.cpp

scoreReads.o:	scoreReads.cu
	nvcc -c scoreReads.cu

clean:
	rm *.o

#scoreGenome.o:  scoreGenome.cpp CScore.h
#	$(CC) -c scoreGenome.cpp

# pattern rule for all objects files
# %.o:  %.cpp %.h
#       $(CC) -c $(input)

# clean:
#       rm -f *~ *.o ; cd ds ; make -f Makefile clean ; cd ..

# cleanall:
#       rm -f *~ *.o *.a; cd ds ; make -f Makefile cleanall ; cd ..
	

