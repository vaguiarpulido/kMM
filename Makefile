CC = g++

HFILES =	CModel.h \
	CSequences.h \
	CScore.h 

CFILES =	CModel.cpp \
	CSequences.cpp 

OFILES =	CModel.o \
	CSequences.o 

test5:	test5.o $(OFILES)
	$(CC) -o test5 $(OFILES) test5.o

test5.o:	CModel.h CSequences.h test5.cpp
	$(CC) -c test5.cpp

CModel.o:	CModel.h CModel.cpp
	$(CC) -c CModel.cpp

CSequences.o:	CSequences.h CSequences.cpp
	$(CC) -c CSequences.cpp

# scoreGenome.o:	scoreGenome.cpp CScore.h
# 	$(CC) -c scoreGenome.cpp

# pattern rule for all objects files
# %.o:	%.cpp %.h
# 	$(CC) -c $(input)

# clean:
# 	rm -f *~ *.o ; cd ds ; make -f Makefile clean ; cd ..

# cleanall:
# 	rm -f *~ *.o *.a; cd ds ; make -f Makefile cleanall ; cd ..


