!--------------------------------------------
!EGraSP: Evolucion GRAvitacional de Sistemas de Particulas
!Maestria en Astrofisica, Universidad de Costa Rica - Jeudy Blanco
!Ultimo cambio: 27/02/2012 -
!--------------------------------------------
program principal
	use Constantes
	use Funciones
	use Auxiliar
	use Tipos
	use Vectores
	use Octree
	use Dinamica
	use Fisica
	use Fuerzas
	use Almacenamiento

	implicit none
	
	include 'mpif.h'

	integer N, Ni, Nproc, i, itr_procs, j, k, i_proc, id_particula_colision, initial_i, save_at, errcode, myid, numprocs , tag, itr_inicio, itr_final, proc_i, n_vecinos, ii, iii
	
	real*8 , allocatable, dimension(:) :: masas, pos_x, pos_y, pos_z, distancias, v_x, v_y, v_z, acc_x, acc_y, acc_z, densidades, densidades_locales, densidades_locales_proc, presiones, presiones_proc
	real*8 , allocatable, dimension(:) :: masas_proc, pos_x_proc, pos_y_proc, pos_z_proc, distancias_proc, v_x_proc, v_y_proc, v_z_proc, acc_x_proc, acc_y_proc, acc_z_proc, densidades_proc
	integer, allocatable, dimension(:,:) :: matriz_vecinos, matriz_vecinos_proc
	integer ipar(0:1)
	real*8 dpar(0:1)
	type(OctreeNode), POINTER :: Arbol
	type(OctreeNode), POINTER ::  NodosParticulas(:)
	real*8 dt, umbralBH, tolerancia_colision, dist_max, beta, temperatura
	type(Particula) p	
	character(len=100) :: path
	logical hubo_colisiones
	character(len=10) :: filename
	character(256) :: namelistfile, prgname
	
	character(len=15) :: tiempoI, tiempoF
	
	character(len=8) :: fechaI, fechaF	
	
	real*8 acc_vect(0:2), posicion_vect(0:2), vector_posicion(0:2)

	integer status(MPI_STATUS_SIZE) 

	integer response(1)
	integer ipunit
	
	namelist /simparam/ N,path,temperatura,umbralBH,initial_i,beta, n_vecinos,dt,j,save_at,tolerancia_colision

	call MPI_INIT(errcode)	
	
	call MPI_COMM_RANK( MPI_COMM_WORLD, myid, errcode )
	call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, errcode )	
	
	write(*,*) "MPI inicializado", " Process ", myid, " of ", numprocs, " is alive"
	!Solamente el root va a pedir los parámetros
	
	!n_vecinos = CANTIDAD_VECINOS
	
	if(myid == 0) then
	
		call getarg(0, prgname)
		call getarg(1, namelistfile)

		ipunit = 100
		open(ipunit, file=namelistfile, status='old', &
			  action='read', err=100)
		read(ipunit, simparam, err=104)
		close(ipunit)

		print *, N,path,temperatura,umbralBH,initial_i,beta, n_vecinos,dt,j,save_at,tolerancia_colision
	endif

	!Transmito los parámetros a los demás nodos
	call MPI_BCAST(N,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(initial_i,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(j,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(save_at,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(umbralBH,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(tolerancia_colision,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(beta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(temperatura,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,errcode)
	call MPI_BCAST(n_vecinos,1,MPI_INTEGER,0,MPI_COMM_WORLD,errcode)

	Ni = N / numprocs

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
	ALLOCATE(distancias(0:N-1))
	ALLOCATE(densidades(0:N-1))		

	ALLOCATE(masas_proc(0:N-1))
	ALLOCATE(pos_x_proc(0:N-1))
	ALLOCATE(pos_y_proc(0:N-1))
	ALLOCATE(pos_z_proc(0:N-1))
	ALLOCATE(v_x_proc(0:N-1))
	ALLOCATE(v_y_proc(0:N-1))
	ALLOCATE(v_z_proc(0:N-1))
	ALLOCATE(acc_x_proc(0:N-1))
	ALLOCATE(acc_y_proc(0:N-1))
	ALLOCATE(acc_z_proc(0:N-1))	
	ALLOCATE(distancias_proc(0:N-1))
	ALLOCATE(densidades_proc(0:N-1))
	ALLOCATE(densidades_locales(0:N-1))
	ALLOCATE(densidades_locales_proc(0:N-1))
	ALLOCATE(presiones(0:N-1))
	ALLOCATE(presiones_proc(0:N-1))
	ALLOCATE(matriz_vecinos(0:N-1,0:n_vecinos-1))
	ALLOCATE(matriz_vecinos_proc(0:N-1,0:n_vecinos-1))	

	ALLOCATE(NodosParticulas(0:N-1))

	!Inicializo la raiz octree
	ALLOCATE(Arbol)

	

	if(myid == 0) then

		!Solo el root va a cargar la nube desde file
		call CargarNube("./data/" // path, masas, pos_x, pos_y, pos_z, v_x, v_y, v_z, distancias, densidades, N)

		dist_max = maxval(distancias)		

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
		Arbol%radio = dist_max		

		call CrearOctree(masas, pos_x, pos_y, pos_z, densidades, N, Arbol, NodosParticulas)		

		!write(*,*) "Octree creado 1era vez (root)"
		
		acc_x = 0.0D+0
		acc_y = 0.0D+0
		acc_z = 0.0D+0
		
		write(*,*) "Recuerde calcular iteraciones y dt basado en freefall time y densidad que hay que calcular"
	endif	
			
	if(myid == 0) then
		tag = 0	
		write(*,*) "Inicia envio"	
		do i = 1, numprocs - 1, 1
			call MPI_SSEND(masas, N, MPI_DOUBLE_PRECISION, i, tag, MPI_COMM_WORLD, errcode)
			call MPI_SSEND(pos_x, N, MPI_DOUBLE_PRECISION, i, tag, MPI_COMM_WORLD, errcode)
			call MPI_SSEND(pos_y, N, MPI_DOUBLE_PRECISION, i, tag, MPI_COMM_WORLD, errcode)
			call MPI_SSEND(pos_z, N, MPI_DOUBLE_PRECISION, i, tag, MPI_COMM_WORLD, errcode)
			call MPI_SSEND(distancias, N, MPI_DOUBLE_PRECISION, i, tag, MPI_COMM_WORLD, errcode)
			call MPI_SSEND(v_x, N, MPI_DOUBLE_PRECISION, i, tag, MPI_COMM_WORLD, errcode)
			call MPI_SSEND(v_y, N, MPI_DOUBLE_PRECISION, i, tag, MPI_COMM_WORLD, errcode)
			call MPI_SSEND(v_z, N, MPI_DOUBLE_PRECISION, i, tag, MPI_COMM_WORLD, errcode)
			call MPI_SSEND(acc_x, N, MPI_DOUBLE_PRECISION, i, tag, MPI_COMM_WORLD, errcode)
			call MPI_SSEND(acc_y, N, MPI_DOUBLE_PRECISION, i, tag, MPI_COMM_WORLD, errcode)
			call MPI_SSEND(acc_z, N, MPI_DOUBLE_PRECISION, i, tag, MPI_COMM_WORLD, errcode)
			call MPI_SSEND(densidades, N, MPI_DOUBLE_PRECISION, i, tag, MPI_COMM_WORLD, errcode)			
		enddo
		write(*,*) "Termina envio"
	else
		write(*,*) "Recibiendo datos desde nodo: ", myid, " N = ", N
		tag = MPI_ANY_TAG
		call MPI_RECV(masas, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
		call MPI_RECV(pos_x, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
		call MPI_RECV(pos_y, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
		call MPI_RECV(pos_z, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
		call MPI_RECV(distancias, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
		call MPI_RECV(v_x, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
		call MPI_RECV(v_y, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
		call MPI_RECV(v_z, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
		call MPI_RECV(acc_x, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
		call MPI_RECV(acc_y, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
		call MPI_RECV(acc_z, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
		call MPI_RECV(densidades, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
		write(*,*) myid, " Done receiving"
	endif
	
	if(myid == 0) then
		tag = MPI_ANY_TAG
		do i = 1, numprocs - 1, 1
			call MPI_RECV(response, 1, MPI_INTEGER, i, tag, MPI_COMM_WORLD, status, errcode)
		enddo
		!write(*,*) "Recibidos todos los mensajes de hijos..."
	else
		response(1) = 0
		tag = 0
		call MPI_SSEND(response, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, errcode)
	endif

	!En este punto, todos los hijos tienen los datos

	!write(*,*) "Continua .... ", myid

	if(.not. myid == 0) then

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

		call CrearOctree(masas, pos_x, pos_y, pos_z, densidades, N, Arbol, NodosParticulas)		

		!write(*,*) "Octree creado 1era vez (hijo : ", myid, ")"

	endif

	!Ciclo principal

	do i = initial_i, j + initial_i, 1
		itr_inicio = myid * Ni
		itr_final = itr_inicio + Ni - 1
		Nproc = Ni
		
		if(numprocs * Ni < N .and.  myid == numprocs - 1) then
			itr_final = itr_final + (N - (numprocs * Ni))			
			Nproc = Nproc + (N - (numprocs * Ni))
		endif
		
		!write(*,*) "Probando rango de particulas a procesar: ", itr_inicio, itr_final, " proceso: ", myid
		
		
		
		call date_and_time (time=tiempoI, date=fechaI)
		call calcularVecinos_Densidad_Presion(N, itr_inicio, itr_final, Arbol, NodosParticulas, masas, pos_x, pos_y, pos_z, v_x, v_y, v_z, n_vecinos, densidades_locales, matriz_vecinos, presiones, temperatura, myid)
		
		!!write(*,*) myid, " - Estoy luego de llamada a calcularVecinos_Densidad_Presion. itr_final = ", itr_final, " - vecinos: ", matriz_vecinos(itr_final, 0:n_vecinos-1)
		
		presiones_proc = presiones(itr_inicio:itr_final)
		densidades_locales_proc = densidades_locales(itr_inicio:itr_final)

		!write(*,*) myid, " Calculadas presiones, densidades y vecinos. Rango: ",  itr_inicio, itr_final, " Presiones: ", presiones_proc

		!---------------------------------------------------
		!---------------------------------------------------				
		!Debo sincronizar presiones y densidades locales antes de continuar con el calculo de aceleracion

		!Recibo presiones locales desde hijos y unifico en root
		!No continua hasta unificarlas
		if(myid == 0) then
			tag = MPI_ANY_TAG
			do itr_procs = 1, numprocs - 1, 1
				
				itr_inicio = itr_procs * Ni
				itr_final = itr_inicio + Ni - 1
				Nproc = Ni
				
				if(numprocs * Ni < N .and.  itr_procs == numprocs - 1) then
					itr_final = itr_final + (N - (numprocs * Ni))
					Nproc = Nproc + (N - (numprocs * Ni))
				endif				
				
				presiones_proc = 0.0D+0
				densidades_locales_proc = 0.0D+0
				
				call MPI_RECV(presiones_proc, Nproc, MPI_DOUBLE_PRECISION, itr_procs, tag, MPI_COMM_WORLD, status, errcode)
				call MPI_RECV(densidades_locales_proc, Nproc, MPI_DOUBLE_PRECISION, itr_procs, tag, MPI_COMM_WORLD, status, errcode)
				
				presiones(itr_inicio:itr_final) = presiones_proc(0:Nproc-1) 
				densidades_locales(itr_inicio:itr_final) = densidades_locales_proc(0:Nproc-1) 
								
			enddo
		else
			!Cada hijo envia las presiones calculadas a padre
			response(1) = 0
			tag = 0		
			
!			write(*,*) "El valor de Nproc es: ", Nproc, " Itr inicio-final: ", itr_inicio, itr_final
			
			presiones_proc = 0.0D+0
			presiones_proc(0:Nproc-1) = presiones(itr_inicio:itr_final)

			densidades_locales_proc = 0.0D+0
			densidades_locales_proc(0:Nproc-1) = densidades_locales(itr_inicio:itr_final)

			
			!write(*,*) myid, " Enviando presiones desde hijo a root: ", presiones_proc
			
			call MPI_SSEND(presiones_proc, Nproc, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, errcode)
			call MPI_SSEND(densidades_locales_proc, Nproc, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, errcode)
		endif	
		
		if(myid == 0) then
			!Envia presiones unificadas a procesos nodos
			tag = 0
			do itr_procs = 1, numprocs - 1, 1
				call MPI_SSEND(presiones, N, MPI_DOUBLE_PRECISION, itr_procs, 0, MPI_COMM_WORLD, errcode)
				call MPI_SSEND(densidades_locales, N, MPI_DOUBLE_PRECISION, itr_procs, 0, MPI_COMM_WORLD, errcode)
			enddo
			!write(*,*) "Enviando presiones a nodos desde root - END - ", presiones
		else
			!write(*,*) "Recibiendo datos desde nodo: ", myid, " N = ", N
			tag = MPI_ANY_TAG
			call MPI_RECV(presiones, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
			call MPI_RECV(densidades_locales, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
		endif			
		
		!---------------------------------------------------
		!---------------------------------------------------
		
		
		itr_inicio = myid * Ni
		itr_final = itr_inicio + Ni - 1
		Nproc = Ni
		
		if(numprocs * Ni < N .and.  myid == numprocs - 1) then
			itr_final = itr_final + (N - (numprocs * Ni))			
			Nproc = Nproc + (N - (numprocs * Ni))
		endif		
		
		!En este punto, todas las presiones deberian estar calculadas
		
		!write(*,*) "-----------------------"
		!write(*,*) myid, " Presiones calculadas?? : ", presiones, " Densidades locales: ", densidades_locales
		!write(*,*) "-----------------------"
		
		!!write(*,*) myid, " - Estoy antes de llamada a pasoLeapFrog. itr_final = ", itr_final, " - vecinos: ", matriz_vecinos(itr_final, 0:n_vecinos-1)
		call pasoLeapFrog(N, itr_inicio, itr_final, Arbol, NodosParticulas, masas, pos_x, pos_y, pos_z, densidades, densidades_locales, v_x, v_y, v_z, acc_x, acc_y, acc_z, dt, umbralBH, tolerancia_colision, beta, n_vecinos, matriz_vecinos, presiones, temperatura, myid)

		if(myid == 0) then
			tag = MPI_ANY_TAG
			do itr_procs = 1, numprocs - 1, 1
				
				itr_inicio = itr_procs * Ni
				itr_final = itr_inicio + Ni - 1
				Nproc = Ni
				
				if(numprocs * Ni < N .and.  itr_procs == numprocs - 1) then
					itr_final = itr_final + (N - (numprocs * Ni))
					Nproc = Nproc + (N - (numprocs * Ni))
				endif				
				
				masas_proc = 0.0D+0
				pos_x_proc = 0.0D+0
				pos_y_proc = 0.0D+0
				pos_z_proc = 0.0D+0
				distancias_proc = 0.0D+0
				v_x_proc = 0.0D+0
				v_y_proc = 0.0D+0
				v_z_proc = 0.0D+0
				acc_x_proc = 0.0D+0
				acc_y_proc = 0.0D+0
				acc_z_proc = 0.0D+0
				densidades_proc = 0.0D+0
				
				call MPI_RECV(masas_proc, Nproc, MPI_DOUBLE_PRECISION, itr_procs, tag, MPI_COMM_WORLD, status, errcode)
				call MPI_RECV(pos_x_proc, Nproc, MPI_DOUBLE_PRECISION, itr_procs, tag, MPI_COMM_WORLD, status, errcode)
				call MPI_RECV(pos_y_proc, Nproc, MPI_DOUBLE_PRECISION, itr_procs, tag, MPI_COMM_WORLD, status, errcode)
				call MPI_RECV(pos_z_proc, Nproc, MPI_DOUBLE_PRECISION, itr_procs, tag, MPI_COMM_WORLD, status, errcode)
				call MPI_RECV(v_x_proc, Nproc, MPI_DOUBLE_PRECISION, itr_procs, tag, MPI_COMM_WORLD, status, errcode)
				call MPI_RECV(v_y_proc, Nproc, MPI_DOUBLE_PRECISION, itr_procs, tag, MPI_COMM_WORLD, status, errcode)
				call MPI_RECV(v_z_proc, Nproc, MPI_DOUBLE_PRECISION, itr_procs, tag, MPI_COMM_WORLD, status, errcode)
				
				!write(*,*) "Recibiendo resultado leapfrog de hijo: ", itr_procs, itr_inicio, itr_final, Nproc
				
				masas(itr_inicio:itr_final) = masas_proc(0:Nproc-1) 
				pos_x(itr_inicio:itr_final) = pos_x_proc(0:Nproc-1) 
				pos_y(itr_inicio:itr_final) = pos_y_proc(0:Nproc-1) 
				pos_z(itr_inicio:itr_final) = pos_z_proc(0:Nproc-1) 
				v_x(itr_inicio:itr_final) = v_x_proc(0:Nproc-1) 
				v_y(itr_inicio:itr_final) = v_y_proc(0:Nproc-1) 
				v_z(itr_inicio:itr_final) = v_z_proc(0:Nproc-1) 
			enddo
			!write(*,*) "Recibidos valores de hijos tras leapfrog..."
		else
			response(1) = 0
			tag = 0
			!call MPI_SSEND(response, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, errcode)
						
			masas_proc(0:Nproc-1) = masas(itr_inicio:itr_final)
			pos_x_proc(0:Nproc-1) = pos_x(itr_inicio:itr_final)
			pos_y_proc(0:Nproc-1) = pos_y(itr_inicio:itr_final)
			pos_z_proc(0:Nproc-1) = pos_z(itr_inicio:itr_final)
			distancias_proc(0:Nproc-1) = distancias(itr_inicio:itr_final)
			v_x_proc(0:Nproc-1) = v_x(itr_inicio:itr_final)
			v_y_proc(0:Nproc-1) = v_y(itr_inicio:itr_final)
			v_z_proc(0:Nproc-1) = v_z(itr_inicio:itr_final)
			acc_x_proc(0:Nproc-1) = acc_x(itr_inicio:itr_final)
			acc_y_proc(0:Nproc-1) = acc_y(itr_inicio:itr_final)
			acc_z_proc(0:Nproc-1) = acc_z(itr_inicio:itr_final)
			densidades_proc(0:Nproc-1) = densidades(itr_inicio:itr_final)			
			tag = 0
			call MPI_SSEND(masas_proc, Nproc, MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_WORLD, errcode)
			call MPI_SSEND(pos_x_proc, Nproc, MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_WORLD, errcode)
			call MPI_SSEND(pos_y_proc, Nproc, MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_WORLD, errcode)
			call MPI_SSEND(pos_z_proc, Nproc, MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_WORLD, errcode)
			call MPI_SSEND(v_x_proc, Nproc, MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_WORLD, errcode)
			call MPI_SSEND(v_y_proc, Nproc, MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_WORLD, errcode)
			call MPI_SSEND(v_z_proc, Nproc, MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_WORLD, errcode)
			
		endif	
				
		!Vuelve a sincronizar con todos los nodos
		
		if(myid == 0) then
			tag = 0
			
			!07/02/2012 Ajusta densidades de particulas a la densidad local para la siguiente iteración
!			densidades(0:N-1) = densidades_locales(0:N-1)

			!write(*,*) "----------------"			
			!write(*,*) "Densidades locales: ", densidades_locales			
			!write(*,*) "Densidades: ", densidades
			!write(*,*) "----------------"
			
			do itr_procs = 1, numprocs - 1, 1
				call MPI_SSEND(masas, N, MPI_DOUBLE_PRECISION, itr_procs, tag, MPI_COMM_WORLD, errcode)
!				call MPI_SSEND(densidades, N, MPI_DOUBLE_PRECISION, itr_procs, tag, MPI_COMM_WORLD, errcode)
				call MPI_SSEND(pos_x, N, MPI_DOUBLE_PRECISION, itr_procs, tag, MPI_COMM_WORLD, errcode)
				call MPI_SSEND(pos_y, N, MPI_DOUBLE_PRECISION, itr_procs, tag, MPI_COMM_WORLD, errcode)
				call MPI_SSEND(pos_z, N, MPI_DOUBLE_PRECISION, itr_procs, tag, MPI_COMM_WORLD, errcode)
				call MPI_SSEND(v_x, N, MPI_DOUBLE_PRECISION, itr_procs, tag, MPI_COMM_WORLD, errcode)
				call MPI_SSEND(v_y, N, MPI_DOUBLE_PRECISION, itr_procs, tag, MPI_COMM_WORLD, errcode)
				call MPI_SSEND(v_z, N, MPI_DOUBLE_PRECISION, itr_procs, tag, MPI_COMM_WORLD, errcode)
			enddo
			!write(*,*) "Enviando datos a nodos desde root - END"
		else
			!write(*,*) "Recibiendo datos desde nodo: ", myid, " N = ", N
			tag = MPI_ANY_TAG
			call MPI_RECV(masas, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
!			call MPI_RECV(densidades, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
			call MPI_RECV(pos_x, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
			call MPI_RECV(pos_y, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
			call MPI_RECV(pos_z, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
			call MPI_RECV(v_x, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
			call MPI_RECV(v_y, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
			call MPI_RECV(v_z, N, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, errcode)
		endif		
		
		call date_and_time (time=tiempoF, date=fechaF)
		
		!En este punto todos los hijos terminaron su parte
		
		!if(myid == 0) then
		if(myid == 0 .and. MOD(i, 10) == 0) then
			write(*,*) "Iteracion ", i, " tiempo transcurrido: ", (i + 1) * dt , "  yrs. Inicio: ", fechaI, " - ", tiempoI, " - Termina: ", fechaF, " - ", tiempoF
		endif

		if(myid == 0 .and. i > 0 .and. MOD(i, save_at) == 0) then									
			Write(filename, '(i10)' )  i
			!write(*,*) "Iteracion ", i, " tiempo transcurrido: ", (i + 1) * dt , "  yrs"
			path = "./data/resultados/" // TRIM(adjustl(filename)) // "_step.csv"						 
			call guardarNube(UNIT_NUBE, path, N, masas, pos_x, pos_y, pos_z, v_x, v_y, v_z, densidades)			
		endif			

	enddo

	if(myid == 0) then
		path = "./data/resultados/complete.csv"			 
		call guardarNube(UNIT_NUBE, path, N, masas, pos_x, pos_y, pos_z, v_x, v_y, v_z, densidades)					
	endif
	
	if(myid == 0) then
		tag = MPI_ANY_TAG
		do i = 1, numprocs - 1, 1
			call MPI_RECV(response, 1, MPI_INTEGER, i, tag, MPI_COMM_WORLD, status, errcode)
		enddo
		!write(*,*) "Recibidos todos los mensajes de hijos..."
	else
		response(1) = 0
		tag = 0
		call MPI_SSEND(response, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, errcode)
	endif	
	
	call MPI_FINALIZE (errcode)
	
	write(*,*) "FIN ", myid

100	write ( 6, * ) 'Cannot read namelist stanza: simparam',  &
		  & trim(namelistfile)
	call MPI_abort (errcode)
104	write ( 6, * ) 'Cannot read namelist stanza: simparam',  &
		  & trim(namelistfile)
	call MPI_abort (errcode)

end
