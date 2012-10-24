program principal
	
	use Almacenamiento
	use Constantes
	use Vectores
	use Funciones
	use Fisica
	use Tipos 
	use BarnesHut
	use Fuerzas
	use Octree

	implicit none
	
	real*8 , allocatable, dimension(:) :: masas, pos_x, pos_y, pos_z, distancias, v_x, v_y, v_z, densidades, distancias_vecinos
	character(len=256) :: path
	character(len=1) :: detalle
	character(len=10) :: tiempoI, tiempoF
	character(len=8) :: fechaI, fechaF
	real*8 ag_bh, ag_fb
	real*8 acc_vect(0:2), acc_vect_total(0:2)
	type(Particula) p
	type(OctreeNode), POINTER :: Arbol, Dummy, Hijo
	type(OctreeNode), POINTER ::  NodosParticulas(:)
	real*8 maxdist, theta
	logical verbose, imprimir
	integer i, j, N, tipo_prueba, n_vecinos, detectados

	integer , allocatable, dimension(:) :: lista_vecinos

	ALLOCATE(Arbol)

!	write(*,*) "Archivo: "	
!	read(*,*) path	

	path = "nubeprueba.csv"
	N = 20

!	write(*,*) "Número de partículas: "	
!	read(*,*) N	

	ALLOCATE(masas(0:N-1))
	ALLOCATE(pos_x(0:N-1))
	ALLOCATE(pos_y(0:N-1))
	ALLOCATE(pos_z(0:N-1))
	ALLOCATE(v_x(0:N-1))
	ALLOCATE(v_y(0:N-1))
	ALLOCATE(v_z(0:N-1))
	ALLOCATE(distancias(0:N-1))
	ALLOCATE(densidades(0:N-1))		
	ALLOCATE(NodosParticulas(0:N-1))

	call CargarNube("./datos/" // path, masas, pos_x, pos_y, pos_z, v_x, v_y, v_z, distancias, densidades, N)

	write(*,*) "Nube cargada"

	maxdist = maxval(distancias)
	Arbol%id = 0
	Arbol%id_particula = -1
	Arbol%hijos_creados = .false.
	Arbol%n_particulas = 0
	Arbol%hoja = .true.
	Arbol%centroide(0) = 0.0D+0
	Arbol%centroide(1) = 0.0D+0
	Arbol%centroide(2) = 0.0D+0
	Arbol%centro_masa(0) = 0.0D+0
	Arbol%centro_masa(1) = 0.0D+0
	Arbol%centro_masa(2) = 0.0D+0
	Arbol%radio = maxdist

	ALLOCATE(Dummy)

	allocate (Dummy%hijo0)
	allocate (Dummy%hijo1)
	allocate (Dummy%hijo2)
	allocate (Dummy%hijo3)
	allocate (Dummy%hijo4)
	allocate (Dummy%hijo5)
	allocate (Dummy%hijo6)
	allocate (Dummy%hijo7)
	Dummy%hijos_creados = .true.

	write(*,*) "Allocated"

	Hijo => DarHijo(Arbol, 1)
	write(*,*) "Hijo cero: ", Hijo%id	
	deallocate(Hijo) 


!	do i = 0, 7, 1
!		Hijo => Dummy%hijos(i)
!		deallocate(Hijo) 
!	enddo	

	deallocate(Dummy)	
	
	write(*,*) "Dellocated"

	stop	

	call CrearOctree(masas, pos_x, pos_y, pos_z, densidades, N, Arbol, NodosParticulas)

	write(*,*) 'Arbol creado'
	
	call imprimirArbol(Arbol)
		write(*,*) "--------"
		write(*,*) "--------"	
		write(*,*) "--------"
		write(*,*) "--------"	
		write(*,*) "--------"
		write(*,*) "--------"	

	call limpiarArbol(Arbol)

	write(*,*) 'Memoria limpia'

	stop
end
