f90 = gfortran
TARGET = Brusselator1D
VPATH = /nfshome/tab7v/QualifyingExam/TestFunctions
MYPATH = /nfshome/tab7v/QualifyingExam/Brusselator1D
OBJECTS = MTT1D.o\
	  inverse.o\
	  IPC.o\

$(TARGET): Brusselator1D.o $(OBJECTS)
	$(f90) -o $(TARGET) $(TARGET).o $(OBJECTS) -L\usr\lib -llapack

$(TARGET).o: Brusselator1D.f90
	$(f90) -c Brusselator1D.f90   

inverse.o: inverse.f90
	$(f90) -I$(VPATH) -c $(VPATH)/inverse.f90 -o $(MYPATH)/inverse.o

MTT1D.o: MTT1D.f90
	$(f90) -I$(VPATH) -c $(VPATH)/MTT1D.f90 -o $(MYPATH)/MTT1D.o

IPC.o: IPC.f90
	$(f90) -I$(VPATH) -c $(VPATH)/IPC.f90 -o $(MYPATH)/IPC.o

.PHONY: clean

clean:
	rm -f $(TARGET)  *.o *~ *.txt
