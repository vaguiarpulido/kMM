CC = g++

HFILES =        CKmers.h \
        CModel.h \
        CSequences.h \
        CScore.h

CFILES =        CKmers.cpp \
        CModel.cpp \
        CSequences.cpp \
	CScore.cu \
        #CScore.cpp \
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

#extractKmers.o:       extractKmers.cpp
#	$(CC) -c extractKmers.cpp

scoreGenome.o:       scoreGenome.cpp
	$(CC) -c scoreGenome.cpp

scoreReads.o:	scoreReads.cu
	nvcc -c scoreReads.cu

CScore.o:	CScore.h CScore.cu
	nvcc -c CScore.cu

#CScore.o:       CScore.h CScore.cpp
#	$(CC) -c CScore.cpp

test5:	test5.o $(OFILES)
	$(CC) -o test5 $(OFILES) test5.o

test5.o:	CModel.h CSequences.h test5.cpp
	$(CC) -c test5.cpp

#CModel.o:	CModel.h CModel.cpp
#	$(CC) -c CModel.cpp

CSequences.o:	CSequences.h CSequences.cpp
	$(CC) -c CSequences.cpp

# scoreGenome.o:	scoreGenome.cpp CScore.h
# 	$(CC) -c scoreGenome.cpp

# pattern rule for all objects files
# %.o:	%.cpp %.h
# 	$(CC) -c $(input)

clean:
	rm -f *~ *.o

# cleanall:
# 	rm -f *~ *.o *.a; cd ds ; make -f Makefile cleanall ; cd ..


