FC = mpif90

FCFLAGS = -ffree-line-length-500 -O2

PROGRAMS = egrasp generar_nube

all: $(PROGRAMS)
	
	
egrasp: almacenamiento.mod fuerzas.mod fisica.mod dinamica.mod octree.mod tipos.mod vectores.mod auxiliar.mod constantes.mod funciones.mod egrasp_ncio.mod
	$(FC) ./src/egrasp.f95 $(FCFLAGS) -I/usr/include -lnetcdf -lnetcdff -o egrasp *.o
	
generar_nube: constantes.mod auxiliar.mod fisica.mod octree.mod vectores.mod tipos.mod funciones.mod almacenamiento.mod egrasp_ncio.mod
	$(FC) ./src/generar_nube.f95 $(FCFLAGS) -I/usr/include -lnetcdf -lnetcdff -o generar_nube *.o
	rm *.mod *.o

constantes.mod: ./src/constantes.f95
	$(FC) $(FCFLAGS) -c ./src/constantes.f95

egrasp_ncio.mod: ./src/egrasp_ncio.f95 tipos.mod
	$(FC) $(FCFLAGS) -I/usr/include -lnetcdf -lnetcdff -c ./src/egrasp_ncio.f95

almacenamiento.mod: vectores.mod constantes.mod egrasp_ncio.mod
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
	rm -rf datos/*
	mkdir datos/resultados
