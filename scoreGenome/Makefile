CC = g++

HFILES =	scoreGenome.h 

CFILES =	scoreGenome.cpp

OFILES =	scoreGenome.o

test5:	test5.o $(OFILES)
	$(CC) -o test5 $(OFILES) test5.o

test5.o:	scoreGenome.h test5.cpp
	$(CC) -c test5.cpp


#CModel.o:	CModel.h CModel.cpp
#	$(CC) -c CModel.cpp

#CSequences.o:	CSequences.h CSequences.cpp
#	$(CC) -c CSequences.cpp


# pattern rule for all objects files
# %.o:	%.cpp %.h
# 	$(CC) -c $(input)

# clean:
# 	rm -f *~ *.o ; cd ds ; make -f Makefile clean ; cd ..

# cleanall:
# 	rm -f *~ *.o *.a; cd ds ; make -f Makefile cleanall ; cd ..


