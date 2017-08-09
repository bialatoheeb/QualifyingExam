f90 = gfortran
TARGET = Brusselator2D
VPATH = /nfshome/tab7v/QualifyingExam/TestFunctions
MYPATH = /nfshome/tab7v/QualifyingExam/Brusselator2D
OBJECTS = MTT2D.o\
	  inverse.o\
	  meshgrid2.o\
	  IPC.o\

$(TARGET): Brusselator2D.o $(OBJECTS)
	$(f90) -o $(TARGET) $(TARGET).o $(OBJECTS) -L\usr\lib -llapack

$(TARGET).o: Brusselator2D.f90
	$(f90) -c Brusselator2D.f90   

inverse.o: inverse.f90
	$(f90) -I$(VPATH) -c $(VPATH)/inverse.f90 -o $(MYPATH)/inverse.o

meshgrid2.o: meshgrid2.f90
	$(f90) -I$(VPATH) -c $(VPATH)/meshgrid2.f90 -o $(MYPATH)/meshgrid2.o

MTT2D.o: MTT2D.f90
	$(f90) -I$(VPATH) -c $(VPATH)/MTT2D.f90 -o $(MYPATH)/MTT2D.o

IPC.o: IPC.f90
	$(f90) -I$(VPATH) -c $(VPATH)/IPC.f90 -o $(MYPATH)/IPC.o

.PHONY: clean

clean:
	rm -f $(TARGET)  *.o *~ *.txt
