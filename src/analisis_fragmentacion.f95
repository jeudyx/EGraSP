program principal
	use Constantes
	use Funciones
	use Auxiliar
	use Tipos
	use Vectores
	use Octree
	use Fisica
	use Almacenamiento
	
	implicit none

	integer N, i, j, k, n_vecinos, id_padre, detectados, acerca
	real*8 , allocatable, dimension(:) :: masas, pos_x, pos_y, pos_z, distancias, v_x, v_y, v_z, acc_x, acc_y, acc_z, densidades, distancias_vecinos, lista_vecinos_real, densidades_locales, presiones
	
	integer , allocatable, dimension(:) ::  indices_vecinos, lista_vecinos
	integer , allocatable, dimension(:,:) ::  matriz_vecinos
		
	real*8 p1(0:2), p2(0:2), p3(0:2), p4(0:2)
	
	logical , allocatable, dimension(:) :: hermanos
	
	type(OctreeNode), POINTER :: Arbol	

	type(OctreeNode), POINTER ::  NodosParticulas(:)
	
	character(len=15) :: tiempoI, tiempoF
	
	character(len=8) :: fechaI, fechaF
		
	character(len=100) :: path
	
	real*8 distancia, densidad_minima, masa_total, masa_acum, volumen_acum, dummy, radio, area

	real*8 dummy_vect(0:2), fuerza(0:2)

	logical seguir	

	print *, "Numero de particulas :"
	read *, N	

	print *, "Archivo con datos de la nube (formato esperado: X,Y,Z,R,Masa) :"
	read *, path	
		
	print *, "Cantidad de vecinos :"
	read *, n_vecinos	

	ALLOCATE(masas(0:N-1))
	ALLOCATE(pos_x(0:N-1))
	ALLOCATE(pos_y(0:N-1))
	ALLOCATE(pos_z(0:N-1))
	ALLOCATE(v_x(0:N-1))
	ALLOCATE(v_y(0:N-1))
	ALLOCATE(v_z(0:N-1))
	ALLOCATE(acc_x(0:N-1))
	ALLOCATE(acc_y(0:N-1))
	ALLOCATE(acc_z(0:N-1))	
	ALLOCATE(indices_vecinos(0:n_vecinos-1))
	ALLOCATE(lista_vecinos(0:n_vecinos-1))
	ALLOCATE(distancias_vecinos(0:n_vecinos-1))
	
	ALLOCATE(densidades_locales(0:N-1))
	ALLOCATE(presiones(0:N-1))
	ALLOCATE(matriz_vecinos(0:N-1, 0:n_vecinos-1))	
	
	ALLOCATE(distancias(0:N-1))
	ALLOCATE(densidades(0:N-1))
	ALLOCATE(NodosParticulas(0:N-1))
	ALLOCATE(hermanos(0:N-1))
	ALLOCATE(Arbol)
	
	masa_total = 0.0D+0
	masa_acum= 0.0D+0 
	volumen_acum= 0.0D+0
	matriz_vecinos = 0
	densidades_locales = 0.0D+0
	presiones = 0.0D+0
	
	hermanos = .false.
	
	call CargarNube("./data/" // path, masas, pos_x, pos_y, pos_z, v_x, v_y, v_z, distancias, densidades, N)

	print *, "Nube cargada"		

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
	Arbol%radio = maxval(distancias)

	detectados = 0
	
	call CrearOctree(masas, pos_x, pos_y, pos_z, densidades, N, Arbol, NodosParticulas)
	
	write(*,*) "Arbol creado"
	
	!call imprimirOctreeStats(Arbol)
		
	call date_and_time (time=tiempoI, date=fechaI)
	
	lista_vecinos = 0

	!i = 4319

	lista_vecinos = 0
	distancias_vecinos = 0.0D+0
	detectados = 0
		
	if(.false.) then
	
		write(*,*) "Desplegando vecinos"

		if(associated(NodosParticulas(i)%adyacente_superior)) then
			write(*,*) "adyacente_superior: ", NodosParticulas(i)%adyacente_superior%id, NodosParticulas(i)%adyacente_superior%hoja, NodosParticulas(i)%adyacente_superior%id_particula, NodosParticulas(i)%adyacente_superior%n_particulas
		endif

		if(associated(NodosParticulas(i)%adyacente_inferior)) then
			write(*,*) "adyacente_inferior: ", NodosParticulas(i)%adyacente_inferior%id, NodosParticulas(i)%adyacente_inferior%hoja, NodosParticulas(i)%adyacente_inferior%id_particula, NodosParticulas(i)%adyacente_inferior%n_particulas
		endif

		if(associated(NodosParticulas(i)%adyacente_izquierdo)) then
			write(*,*) "adyacente_izquierdo: ", NodosParticulas(i)%adyacente_izquierdo%id, NodosParticulas(i)%adyacente_izquierdo%hoja, NodosParticulas(i)%adyacente_izquierdo%id_particula, NodosParticulas(i)%adyacente_izquierdo%n_particulas
		endif

		if(associated(NodosParticulas(i)%adyacente_derecho)) then
			write(*,*) "adyacente_derecho: ", NodosParticulas(i)%adyacente_derecho%id, NodosParticulas(i)%adyacente_derecho%hoja, NodosParticulas(i)%adyacente_derecho%id_particula, NodosParticulas(i)%adyacente_derecho%n_particulas
		endif

		if(associated(NodosParticulas(i)%adyacente_frente)) then
			write(*,*) "adyacente_frente: ", NodosParticulas(i)%adyacente_frente%id, NodosParticulas(i)%adyacente_frente%hoja, NodosParticulas(i)%adyacente_frente%id_particula, NodosParticulas(i)%adyacente_frente%n_particulas
		endif

		if(associated(NodosParticulas(i)%adyacente_trasero)) then
			write(*,*) "adyacente_trasero: ", NodosParticulas(i)%adyacente_trasero%id, NodosParticulas(i)%adyacente_trasero%hoja, NodosParticulas(i)%adyacente_trasero%id_particula, NodosParticulas(i)%adyacente_trasero%n_particulas
		endif

		write(*,*) "Desplegando hermanos"

		do j = 0, 7, 1
			write(*,*) "Hermano: ", j, NodosParticulas(i)%padre%hijos(j)%id, NodosParticulas(i)%padre%hijos(j)%hoja, NodosParticulas(i)%padre%hijos(j)%id_particula, NodosParticulas(i)%padre%hijos(j)%n_particulas
		enddo

		write(*,*) "-------------------"
	
	endif

	if(.true.) then

		print *, "Digite particula i :"
		read *, i	

		lista_vecinos = 0
		distancias_vecinos = 0.0D+0
		detectados = 0

		call ResetearVisitados(Arbol)
		call Vecinos(NodosParticulas(i), NodosParticulas, NodosParticulas(i), N, detectados, n_vecinos, lista_vecinos, distancias_vecinos, 0, 2, 0)

		do j = 0, n_vecinos - 1, 1		
			indices_vecinos(j) = j
		enddo
		
		dummy = calcularDensidadLocal(N, i, lista_vecinos, n_vecinos, pos_x, pos_y, pos_z, masas)

		write(*,*) "Densidad local de i = ", i, " = ", dummy

		distancias_vecinos = ordenar(distancias_vecinos, indices_vecinos, n_vecinos)

!		do j = 0, n_vecinos - 1, 1		
!			p1(0) = pos_x(i)
!			p1(1) = pos_y(i)
!			p1(2) = pos_z(i)
			
!			p2(0) = pos_x(indices_vecinos(j))
!			p2(1) = pos_y(indices_vecinos(j))
!			p2(2) = pos_z(indices_vecinos(j))

!			p3(0) = v_x(i)
!			p3(1) = v_y(i)
!			p3(2) = v_z(i)

!			p4(0) = v_x(indices_vecinos(j))
!			p4(1) = v_y(indices_vecinos(j))
!			p4(2) = v_z(indices_vecinos(j))
			
!			write(*,*) i, j, " Vab . Rab : ", productoEscalar3D(diferenciaVectores3D(p1, p2), diferenciaVectores3D(p3, p4))
!		enddo

		
		write(*,*) "-------------- VECINOS --------------"

		do j = 0, n_vecinos - 1, 1		
			
			p1(0) = pos_x(i)
			p1(1) = pos_y(i)
			p1(2) = pos_z(i)
			
			p2(0) = pos_x(indices_vecinos(j))
			p2(1) = pos_y(indices_vecinos(j))
			p2(2) = pos_z(indices_vecinos(j))

			p3(0) = v_x(i)
			p3(1) = v_y(i)
			p3(2) = v_z(i)

			p4(0) = v_x(indices_vecinos(j))
			p4(1) = v_y(indices_vecinos(j))
			p4(2) = v_z(indices_vecinos(j))
			
			if(productoEscalar3D(diferenciaVectores3D(p1, p2), diferenciaVectores3D(p3, p4)) < 0) then
				acerca = -1
			else
				acerca = 0
			endif
			
			write(*,*) lista_vecinos(indices_vecinos(j)),",",pos_x(lista_vecinos(indices_vecinos(j))),",",pos_y(lista_vecinos(indices_vecinos(j))),",",pos_z(lista_vecinos(indices_vecinos(j))),",",distancias_vecinos(j), ",", acerca
		enddo

		write(*,*) "-------------Fuerza bruta---------------"

		lista_vecinos = 0

		call VecinosFuerzaBruta(NodosParticulas, i, N, n_vecinos, lista_vecinos)
				
		do j = 0, n_vecinos - 1, 1		
			write(*,*) lista_vecinos(j),",",pos_x(lista_vecinos(j)),",",pos_y(lista_vecinos(j)),",",pos_z(lista_vecinos(j)),",",distanciaPuntos(NodosParticulas(i)%centro_masa, NodosParticulas(lista_vecinos(j))%centro_masa)
		enddo		

		return

	endif
			
	do i = 0, N - 1, 1

		!write(*,*) "i: ", i
	
		lista_vecinos = 0
		distancias_vecinos = 0.0D+0
		detectados = 0
		call ResetearVisitados(Arbol)
		
		do j = 0, N - 1, 1
			NodosParticulas(i)%visitado = .false.
		enddo
		
		call Vecinos(NodosParticulas(i), NodosParticulas, NodosParticulas(i), N, detectados, n_vecinos, lista_vecinos, distancias_vecinos, 0, 2, 0)
		densidades_locales(i) = calcularDensidadLocal(N, i, lista_vecinos, n_vecinos, pos_x, pos_y, pos_z, masas)
		presiones(i) = PresionGasIdeal(densidades_locales(i), 10.0D+0)
		matriz_vecinos(i, 0:n_vecinos - 1) = lista_vecinos

		if(.true. .or. MOD(i, 100) == 0) then

			!write(*,*) "Iteracion: ", i, " densidad local,", densidades_locales(i), " presion: ", presiones(i), presiones(i) / K_BOLTZMANN
			
			!do j = 0, n_vecinos - 1, 1		
			!	indices_vecinos(j) = j
			!enddo

			!distancias_vecinos = ordenar(distancias_vecinos, indices_vecinos, n_vecinos)

			!write(*,*) "id,x,y,z,d"

			!do j = 0, n_vecinos - 1, 1		
			!	write(*,*) lista_vecinos(indices_vecinos(j)),",",pos_x(lista_vecinos(indices_vecinos(j))),",",pos_y(lista_vecinos(indices_vecinos(j))),",",pos_z(lista_vecinos(indices_vecinos(j))),",",distancias_vecinos(j)
			!enddo
		endif	
	enddo
	
	return
	
	do i = 0, N - 1, 1
		lista_vecinos = matriz_vecinos(i, 0:n_vecinos-1)
		radio = calcularRadio(masas(i), densidades(i)) / 100.0D+0	!en metros
		area = 4 * PI * (radio**2)
		dummy_vect = GradientePresion(i, n_vecinos, lista_vecinos, presiones, N, pos_x, pos_y, pos_z, masas, densidades)		
		fuerza = dummy_vect * area				
		write(*,*) "i = ", i, " presion en i: ", presiones(i), 10.0D+0 * presiones(i) / K_BOLTZMANN, " - densidad local: ", densidades_locales(i), " Gradiente P: ", dummy_vect, magnitudVector3D(dummy_vect), " Fuerza grad presion: ", magnitudVector3D(fuerza), " aceleracion: ", fuerza / (masas(i) * SOLAR_MASS_KG), magnitudVector3D(fuerza / (masas(i) * SOLAR_MASS_KG))
	enddo
	
	call date_and_time (time=tiempoF, date=fechaF)

	write(*,*) "Inicio Aproximacion: ", fechaI, " - ", tiempoI, " - Termina Aproximacion: ", fechaF, " - ", tiempoF

	write(*,*) "----------------------"

	return

end
