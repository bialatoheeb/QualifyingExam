f90 = gfortran
TARGET = Brusselator3D
VPATH = /nfshome/tab7v/QualifyingExam/TestFunctions
MYPATH = /nfshome/tab7v/QualifyingExam/Brusselator3D
OBJECTS = MTT3D.o  inverse.o  meshgrid2.o  meshgrid3.o IPC.o

$(TARGET): Brusselator3D.o $(OBJECTS)
	$(f90) -o $(TARGET) $(TARGET).o $(OBJECTS) -L\usr\lib -llapack

$(TARGET).o: Brusselator3D.f90
	$(f90) -c Brusselator3D.f90   

meshgrid2.o: meshgrid2.f90
	$(f90) -I$(VPATH) -c $(VPATH)/meshgrid2.f90 -o $(MYPATH)/meshgrid2.o

meshgrid3.o: meshgrid3.f90
	$(f90) -I$(VPATH) -c $(VPATH)/meshgrid3.f90 -o $(MYPATH)/meshgrid3.o

inverse.o: inverse.f90
	$(f90) -I$(VPATH) -c $(VPATH)/inverse.f90 -o $(MYPATH)/inverse.o

MTT3D.o: MTT3D.f90
	$(f90) -I$(VPATH) -c $(VPATH)/MTT3D.f90 -o $(MYPATH)/MTT3D.o

IPC.o: IPC.f90
	$(f90) -I$(VPATH) -c $(VPATH)/IPC.f90 -o $(MYPATH)/IPC.o

.PHONY: clean

clean:
	rm -f $(TARGET)  *.o *~ *.txt
