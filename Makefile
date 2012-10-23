FC = mpif90

FCFLAGS = -ffree-line-length-1000 -m64 -g

PROGRAMS = prueba_mem egrasp generar_nube

all: $(PROGRAMS)
	
prueba_mem: almacenamiento.mod fuerzas.mod fisica.mod octree.mod constantes.mod funciones.mod
	$(FC) ./src/prueba_mem.f95 $(FCFLAGS)  -o prueba_mem *.o
	
egrasp: almacenamiento.mod fuerzas.mod fisica.mod dinamica.mod octree.mod tipos.mod vectores.mod auxiliar.mod constantes.mod funciones.mod 
	$(FC) ./src/egrasp.f95 $(FCFLAGS)  -o egrasp *.o
	
generar_nube: constantes.mod auxiliar.mod fisica.mod octree.mod vectores.mod tipos.mod funciones.mod almacenamiento.mod 
	$(FC) ./src/generar_nube.f95 $(FCFLAGS)  -o generar_nube *.o
	rm *.mod *.o

constantes.mod: ./src/constantes.f95
	$(FC) $(FCFLAGS) -c ./src/constantes.f95

almacenamiento.mod: vectores.mod constantes.mod 
	$(FC) $(FCFLAGS) -c ./src/almacenamiento.f95

tipos.mod: ./src/tipos.f95
	$(FC) $(FCFLAGS) -c ./src/tipos.f95

auxiliar.mod: tipos.mod
	$(FC) $(FCFLAGS) -c ./src/auxiliar.f95

vectores.mod: tipos.mod constantes.mod
	$(FC) $(FCFLAGS) -c ./src/vectores.f95

funciones.mod: tipos.mod constantes.mod vectores.mod
	$(FC) $(FCFLAGS) -c ./src/funciones.f95

fisica.mod: tipos.mod constantes.mod vectores.mod auxiliar.mod funciones.mod
	$(FC) $(FCFLAGS) -c ./src/fisica.f95

octree.mod: tipos.mod constantes.mod vectores.mod auxiliar.mod funciones.mod fisica.mod
	$(FC) $(FCFLAGS) -c ./src/octree.f95

barneshut.mod: tipos.mod constantes.mod vectores.mod auxiliar.mod funciones.mod fisica.mod octree.mod
	$(FC) $(FCFLAGS) -c ./src/barneshut.f95

fuerzas.mod: barneshut.mod tipos.mod octree.mod
	$(FC) $(FCFLAGS) -c ./src/fuerzas.f95

dinamica.mod: barneshut.mod tipos.mod constantes.mod vectores.mod auxiliar.mod funciones.mod fisica.mod octree.mod fuerzas.mod
	$(FC) $(FCFLAGS) -c ./src/dinamica.f95


#------------------


# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.mod *.o
	rm -f *~ $(PROGRAMS)
#	rm -rf datos/*
#	mkdir datos/resultados
