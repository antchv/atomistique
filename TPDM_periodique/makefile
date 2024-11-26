CXX = g++
CXXFLAGS = -Wall -g -O3 -pg

# Définition de la variable METHODE en fonction de l'option choisie
ifeq ($(METHODE), EULER)
METHODE_FLAG = -DEULER
else ifeq ($(METHODE), VV)
METHODE_FLAG = -DVV
else
METHODE_FLAG = -DVV  # Velocity-Verlet Par défaut
endif

a.out : main.o fcts.o methode.o util.o
	$(CXX) $(CXXFLAGS) $(METHODE_FLAG) main.o fcts.o methode.o util.o -o a.exe

fcts.o : fcts.cpp fcts.h
	$(CXX) $(CXXFLAGS) $(METHODE_FLAG) -c fcts.cpp

methode.o : methode.cpp methode.h
	$(CXX) $(CXXFLAGS) $(METHODE_FLAG) -c methode.cpp

util.o : util.cpp util.h
	$(CXX) $(CXXFLAGS) $(METHODE_FLAG) -c util.cpp

main.o : main.cpp fcts.h methode.h util.h
	$(CXX) $(CXXFLAGS) $(METHODE_FLAG) -c main.cpp

.PHONY: clean

clean:
	/bin/rm -f a.exe *.o
