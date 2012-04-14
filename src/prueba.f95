
logical function imprimir(i)

	implicit none

	integer i

	!if(i == 77 .or. i == 1670 .or. i == 2440 .or. i == 3331 .or. i == 4021 .or. i == 5198 .or. i == 6366 .or. i == 7175 .or. i == 8017 .or. i == 9127) then
	if(i == 77 .or. i == 2440 .or. i == 4021 .or. i == 6366 .or. i == 9127) then
		imprimir = .true.
	else
		imprimir = .false.
	endif

end function

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
	
	!real*8 masas(0:N-1), pos_x(0:N-1), pos_y(0:N-1), pos_z(0:N-1), distancias(0:N-1), v_x(0:N-1), v_y(0:N-1), v_z(0:N-1), densidades(0:N-1)
	real*8 , allocatable, dimension(:) :: masas, pos_x, pos_y, pos_z, distancias, v_x, v_y, v_z, densidades, distancias_vecinos
	character(len=256) :: path
	character(len=1) :: detalle
	character(len=10) :: tiempoI, tiempoF
	character(len=8) :: fechaI, fechaF
	real*8 ag_bh, ag_fb
	real*8 acc_vect(0:2), acc_vect_total(0:2)
	type(Particula) p
	type(OctreeNode), POINTER :: Arbol	
	type(OctreeNode), POINTER ::  NodosParticulas(:)
	real*8 maxdist, theta
	logical verbose, imprimir
	integer i, j, N, tipo_prueba, n_vecinos, detectados

	integer , allocatable, dimension(:) :: lista_vecinos

	ALLOCATE(Arbol)

	tiempoI = '\0'
	tiempoF = '\0'
	fechaI = '\0'
	fechaF = '\0'
	path = '\0'
	acc_vect = 0.0D+0
	acc_vect_total = 0.0D+0

	write(*,*) "Archivo: "	
	read(*,*) path	

	write(*,*) "Número de partículas: "	
	read(*,*) N	

	write(*,*) "Theta:"
	read(*,*), theta

	write(*,*) "Imprimir detalle (y/n): "	
	read(*,*) detalle	

	write(*,*) "Tipo de prueba? 1 = gravitacional, 2 = vecinos"
	read(*,*) tipo_prueba

	if(tipo_prueba == 2) then
		write(*,*) "Cantidad de vecinos a buscar: "	
		read(*,*) n_vecinos	
	endif

	if(detalle == 'y') then
		verbose = .true.
	else
		verbose = .false.
	endif

!	write(*,*) "Parametros: ", path, N, verbose

	!path = "./datos/prueba_performance10k.csv"

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
	ALLOCATE(lista_vecinos(0:n_vecinos-1))
	ALLOCATE(distancias_vecinos(0:n_vecinos-1))

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


	call CrearOctree(masas, pos_x, pos_y, pos_z, densidades, N, Arbol, NodosParticulas)

	write(*,*) "Arbol creado. Distancia máxima: ", maxdist
	
	tiempoI = '\0'
	tiempoF = '\0'
	fechaI = '\0'
	fechaF = '\0'

	if(tipo_prueba == 1) then

		p = construirParticula(0, masas(0), pos_x(0), pos_y(0), pos_z(0), v_x(0), v_y(0), v_z(0), densidades(0))

		call date_and_time (time=tiempoI, date=fechaI)
		acc_vect = calcularAceleracionGFB(p, pos_x, pos_y, pos_z, masas, N)
		call date_and_time (time=tiempoF, date=fechaF)
		ag_fb = magnitudVector3D(acc_vect)
		write(*,*) "Aceleración gravitacional fuerza bruta. Valor: ", acc_vect, ag_fb
		tiempoI = '\0'
		tiempoF = '\0'
		fechaI = '\0'
		fechaF = '\0'
		acc_vect = 0.0D+0
		write(*,*) "------------------------"
		call date_and_time (time=tiempoI, date=fechaI)
		acc_vect = calcularAceleracion(p, Arbol, theta, 0.0D+0, N)
		call date_and_time (time=tiempoF, date=fechaF)
		ag_bh = magnitudVector3D(acc_vect)
		write(*,*) "Aceleración gravitacional barnes-hut. Valor: ", acc_vect, ag_bh
		write(*,*) "------------------------"
		write(*,*) "------------------------"
		call date_and_time (time=tiempoI, date=fechaI)
		do i = 0, N-1, 1
			p = construirParticula(i, masas(i), pos_x(i), pos_y(i), pos_z(i), v_x(i), v_y(i), v_z(i), densidades(i))
			acc_vect = calcularAceleracionGFB(p, pos_x, pos_y, pos_z, masas, N)
			if(verbose .and. imprimir(i)) then		
				write(*,*) i, ",", acc_vect, ",", magnitudVector3D(acc_vect), ","
			endif
	!		acc_vect_total(0) = acc_vect_total(0) + acc_vect(0)
	!		acc_vect_total(1) = acc_vect_total(1) + acc_vect(1)
	!		acc_vect_total(2) = acc_vect_total(2) + acc_vect(2)
		enddo
		call date_and_time (time=tiempoF, date=fechaF)
		write(*,*) "Aceleración gravitacional fuerza bruta. Inicia: ", tiempoI, ", termina: ", tiempoF!, " - Total: ", magnitudVector3D(acc_vect_total)
		write(*,*) "------------------------"
		write(*,*) "------------------------"
		acc_vect_total = 0.0D+0
		call date_and_time (time=tiempoI, date=fechaI)
		do i = 0, N-1, 1
			p = construirParticula(i, masas(i), pos_x(i), pos_y(i), pos_z(i), v_x(i), v_y(i), v_z(i), densidades(i))
			acc_vect = calcularAceleracion(p, Arbol, theta, 0.0D+0, N)
			if(verbose .and. imprimir(i)) then		
				write(*,*) i, ",", acc_vect, ",", magnitudVector3D(acc_vect), ","
			endif
	!		acc_vect_total(0) = acc_vect_total(0) + acc_vect(0)
	!		acc_vect_total(1) = acc_vect_total(1) + acc_vect(1)
	!		acc_vect_total(2) = acc_vect_total(2) + acc_vect(2)
		enddo
		call date_and_time (time=tiempoF, date=fechaF)
		write(*,*) "Aceleración gravitacional Barnes-Hut. Inicia: ", tiempoI, ", termina: ", tiempoF!, " - Total: ", magnitudVector3D(acc_vect_total)
	
	else
!		call date_and_time (time=tiempoI, date=fechaI)
!		do i = 0, N-1, 1
!			lista_vecinos(:) = 0
!			call VecinosFuerzaBruta(NodosParticulas, i, N, n_vecinos, lista_vecinos)
!			if(verbose .and. i == 499) then
!				write(*,*) "Vecinos por fuerza bruta"
!				do j = 0, n_vecinos - 1, 1
!					write(*,*) lista_vecinos(j)
!				enddo
!			endif
!		enddo
!		call date_and_time (time=tiempoF, date=fechaF)
!		write(*,*) "Vecinos por fuerza bruta. Inicia: ", tiempoI, ", termina: ", tiempoF		
		!--------------
		call date_and_time (time=tiempoI, date=fechaI)
		do i = 0, N-1, 1
			detectados = 0
			lista_vecinos(:) = 0
			distancias_vecinos(:) = 0.0D+0
			call ResetearVisitados(Arbol)
		
			do j = 0, N - 1, 1
				NodosParticulas(j)%visitado = .false.
			enddo
			call Vecinos(NodosParticulas(i), NodosParticulas, NodosParticulas(i), N, detectados, n_vecinos, lista_vecinos, distancias_vecinos, 0, 2, 0)
			if(verbose .and. i == 499) then
				write(*,*) "Vecinos por fuerza bruta"
				do j = 0, n_vecinos - 1, 1
					write(*,*) lista_vecinos(j)
				enddo
			endif
		enddo
		call date_and_time (time=tiempoF, date=fechaF)
		write(*,*) "Vecinos por aproximación. Inicia: ", tiempoI, ", termina: ", tiempoF		
	endif

	stop
end
