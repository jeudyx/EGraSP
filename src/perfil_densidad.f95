program principal
	use Constantes
	use Funciones
	use Auxiliar
	use Tipos
	use Vectores
	use Octree
	use Fisica
	use Almacenamiento
	use Fuerzas
	
	implicit none

	integer N, i, j, k, n_vecinos, id_padre, detectados
	real*8 , allocatable, dimension(:) :: masas, pos_x, pos_y, pos_z, distancias, v_x, v_y, v_z, densidades, distancias_vecinos, lista_vecinos_real, densidades_locales, presiones, presiones_netas
	
	integer , allocatable, dimension(:) ::  indices_particulas, lista_vecinos
	integer , allocatable, dimension(:,:) ::  matriz_vecinos
		
	real*8 centro_masa(0:2), gradiente_presion(0:2), acc_vect(0:2), acc_vect_visc(0:2)
	
	real*8 p1(0:2), p2(0:2), p3(0:2), p4(0:2)
		
	type(OctreeNode), POINTER :: Arbol	

	type(OctreeNode), POINTER ::  NodosParticulas(:)
			
	character(len=100) :: path
	
	type(Particula) p
	
	real*8 temperatura, masa_acum, viscosity

	print *, "Numero de particulas :"
	read *, N	

	print *, "Archivo con datos de la nube (formato esperado: X,Y,Z,R,Masa) :"
	read *, path	
		
	print *, "Cantidad de vecinos :"
	read *, n_vecinos	

	print *, "Temperatura :"
	read *, temperatura

	ALLOCATE(masas(0:N-1))
	ALLOCATE(pos_x(0:N-1))
	ALLOCATE(pos_y(0:N-1))
	ALLOCATE(pos_z(0:N-1))
	ALLOCATE(v_x(0:N-1))
	ALLOCATE(v_y(0:N-1))
	ALLOCATE(v_z(0:N-1))
	ALLOCATE(indices_particulas(0:N-1))
	ALLOCATE(lista_vecinos(0:n_vecinos-1))
	ALLOCATE(distancias_vecinos(0:n_vecinos-1))
	
	ALLOCATE(densidades_locales(0:N-1))
	ALLOCATE(presiones(0:N-1))
	ALLOCATE(presiones_netas(0:N-1))	
	ALLOCATE(matriz_vecinos(0:N-1, 0:n_vecinos-1))	
	
	ALLOCATE(distancias(0:N-1))
	ALLOCATE(densidades(0:N-1))
	ALLOCATE(NodosParticulas(0:N-1))
	ALLOCATE(Arbol)
	
	matriz_vecinos = 0
	densidades_locales = 0.0D+0
	presiones = 0.0D+0
		
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
	
	lista_vecinos = 0
	distancias_vecinos = 0.0D+0
	detectados = 0
		
	if(.false.) then

	endif
	
	centro_masa(0) = pos_x(0)
	centro_masa(1) = pos_y(0)
	centro_masa(2) = pos_z(0)
	masa_acum = masas(0)
	
	do i = 1, N - 1, 1		
		centro_masa = centroMasa(masas(i), pos_x(i), pos_y(i), pos_z(i), masa_acum, centro_masa(0), centro_masa(1), centro_masa(2))
		masa_acum = masa_acum + masas(i)
	enddo	
	
	write(*,*) "Centro de masa de conjunto: ", centro_masa
	
	do i = 0, N - 1, 1
		indices_particulas(i) = i
	enddo
	
	distancias = ordenar(distancias, indices_particulas, N)
	
	write(*,*) "i,distancia,distancia_CentroMasa,x,y,z,masa,densidad_local,grad_presion_x,grad_presion_y,grad_presion_z,accp_x,accp_y,accp_z,acc_grav_x,acc_grav_y,acc_grav_z,presion,presion_viscosidad,acc_visc_x,acc_visc_y,acc_visc_z"
	
	do i = 0, N - 1, 1

		!write(*,*) "i: ", i
	
		lista_vecinos = 0
		distancias_vecinos = 0.0D+0
		detectados = 0
		call ResetearVisitados(Arbol)
		
		do j = 0, N - 1, 1
			NodosParticulas(j)%visitado = .false.
		enddo		
				
		call Vecinos(NodosParticulas(indices_particulas(i)), NodosParticulas, NodosParticulas(indices_particulas(i)), N, detectados, n_vecinos, lista_vecinos, distancias_vecinos, 0, 2, 0)
		
		!write(*,*) "Vecinos de i: ", indices_particulas(i), " -- ", lista_vecinos
		
		densidades_locales(indices_particulas(i)) = calcularDensidadLocal(N, indices_particulas(i), lista_vecinos, n_vecinos, pos_x, pos_y, pos_z, masas)
		presiones(indices_particulas(i)) = PresionGasIdeal(densidades_locales(indices_particulas(i)), temperatura)
		matriz_vecinos(indices_particulas(i), 0:n_vecinos - 1) = lista_vecinos
		
		p1(0) = pos_x(indices_particulas(i))
		p1(1) = pos_y(indices_particulas(i))
		p1(2) = pos_z(indices_particulas(i))
		
		!write(*,*) indices_particulas(i), ",", distancias(i), ",", distanciaPuntos(centro_masa, p1), ",", pos_x(indices_particulas(i)), ",", pos_y(indices_particulas(i)), ",", pos_z(indices_particulas(i)), ",", masas(indices_particulas(i)), ",", densidades_locales(indices_particulas(i)), ",", presiones(indices_particulas(i))
	enddo
	
	!write(*,*) "MASAS: ", masas
	!write(*,*) "DENSIDADES: ", densidades
	
	!do i = 1, N - 1, 1		
	!	viscosity = ArtificialViscosity(indices_particulas(i), temperatura, n_vecinos, matriz_vecinos(indices_particulas(i), 0:n_vecinos - 1), N, pos_x, pos_y, pos_z, v_x, v_y, v_z, densidades, densidades_locales, masas)
	!	presiones_netas(indices_particulas(i)) = presiones(indices_particulas(i)) + viscosity
	!enddo		
	
	do i = 0, N - 1, 1
		gradiente_presion = GradientePresion(indices_particulas(i), n_vecinos, matriz_vecinos(indices_particulas(i), 0:n_vecinos - 1), presiones, N, pos_x, pos_y, pos_z, masas, densidades)
		
		acc_vect_visc = ArtificialViscosityAcc(indices_particulas(i), temperatura, n_vecinos, matriz_vecinos(indices_particulas(i), 0:n_vecinos - 1), N, pos_x, pos_y, pos_z, v_x, v_y, v_z, densidades, densidades_locales, masas)
		
		viscosity = ArtificialViscosity(indices_particulas(i), temperatura, n_vecinos, matriz_vecinos(indices_particulas(i), 0:n_vecinos - 1), N, pos_x, pos_y, pos_z, v_x, v_y, v_z, densidades, densidades_locales, masas)
		
		p = construirParticula(indices_particulas(i), masas(indices_particulas(i)), pos_x(indices_particulas(i)), pos_y(indices_particulas(i)), pos_z(indices_particulas(i)), v_x(indices_particulas(i)), v_y(indices_particulas(i)), v_z(indices_particulas(i)), densidades(indices_particulas(i)))
		acc_vect = calcularAceleracion(p, Arbol, 0.5D+0, densidades, N)
		write(*,*) indices_particulas(i), ",", distancias(i), ",", distanciaPuntos(centro_masa, p%posicion), ",", pos_x(indices_particulas(i)), ",", pos_y(indices_particulas(i)), ",", pos_z(indices_particulas(i)), ",", masas(indices_particulas(i)), ",", densidades_locales(indices_particulas(i)), ",", gradiente_presion(0), ",", gradiente_presion(1), ",", gradiente_presion(2), ",", gradiente_presion(0) / (densidades(indices_particulas(i)) / 1000.0D+0), ",", gradiente_presion(1) / (densidades(indices_particulas(i)) / 1000.0D+0), ",", gradiente_presion(2) / (densidades(indices_particulas(i)) / 1000.0D+0), ",", acc_vect(0), ",", acc_vect(1), ",", acc_vect(2), ",", presiones(indices_particulas(i)), ",", viscosity, ",", acc_vect_visc(0), ",", acc_vect_visc(1), ",", acc_vect_visc(2)
	enddo
		
	
	return
end
