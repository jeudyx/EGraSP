!---------------------------------------------------------------------------------------
!---  Modelado de nube molecular, construccion inicial 				     ---|
!---  Maestria en Astrofisica - UCR					             ---|
!---  Jeudy Blanco - A57060                  					     ---|
!---------------------------------------------------------------------------------------

subroutine PerturbarDensidad(N, pos_x, pos_y, pos_z, masas, densidad, n_perturbacion, radio, n_vecinos, n_densidad)
	use Constantes
	use Auxiliar
	use Fisica
	use Octree
	use Vectores
	use Tipos
	use Funciones 
	
	implicit none
	
	integer N, i, idx, itr, j, k, ipart, n_perturbacion, n_vecinos, detectados	
	real*8 pos_x(0:N-1), pos_y(0:N-1), pos_z(0:N-1), masas(0:N-1), distancias(0:N-1), vector_posicion(0:2), vector_posicion_part(0:2), densidades(0:N-1)
	real*8 distancias_vecinos(0:n_vecinos-1), distancias_vecinos_ord(0:n_vecinos-1)
	real*8 radio, radio_particula, densidad_local, densidad, randomX, randomY, randomZ, segmento, inicio, final, posicion, distancia, n_densidad
	real*8 nueva_posicion(0:2)
	integer indices(0:N-1), lista_vecinos(0:n_vecinos-1), indices_vecinos(0:n_vecinos-1)
	logical visitados(0:N-1), condicion_densidad, cambio_signo, es_vecino
	type(OctreeNode), POINTER :: Arbol	
	type(OctreeNode) NodosParticulas(0:N-1)

	ALLOCATE(Arbol)
	
	segmento = (radio*2.0D+0) / n_perturbacion
	
	write(*,*) "Tamaño de segmento: ", segmento
	
	cambio_signo = .false.
		
	densidades = densidad
			
	visitados = .false.			
		
	ipart = 1

	inicio = -1.0D+0 * radio
	final = inicio + segmento
	
	posicion = (inicio + final) / 2.0D+0
	
	write(*,*) "INICIA: ", inicio, final, posicion
	
	do idx = 0, n_perturbacion - 1, 1
		
		pos_x(idx) = posicion
		
		inicio = final
		final = final + segmento
		posicion = (inicio + final) / 2.0D+0
		
		pos_y(idx) = 0.0D+0
		pos_z(idx) = 0.0D+0
	
		!Tiene que reacomodar partículas cerca a particula en puntos de perturbacion hasta que densidad local sea 10 veces la promedio
		
		CALL init_random_seed()
		
		condicion_densidad = .false.
				
		visitados(idx) = .true.
		
		!radio_particula = calcularRadio(masas(idx), densidad) / PARSEC_CMS
				
		!write(*,*) "HOLA"
		
		write(*,*) "Posicion particula en zona de perturbacion: ", pos_x(idx), pos_y(idx), pos_z(idx)
		
		condicion_densidad = .false.
	enddo
		
	do idx = 0, n_perturbacion - 1, 1
		
		vector_posicion_part(0) = pos_x(idx)
		vector_posicion_part(1) = pos_y(idx)
		vector_posicion_part(2) = pos_z(idx)
		
		CALL init_random_seed()
		
		condicion_densidad = .false.

		do i = 0, N - 1, 1
			vector_posicion(0) = pos_x(i)
			vector_posicion(1) = pos_y(i)
			vector_posicion(2) = pos_z(i)

			distancias(i) = magnitudVector3D(vector_posicion)
			indices(i) = i
		enddo	

		call limpiarArbol(Arbol)
		Arbol%radio = maxval(distancias)
		call CrearOctree(masas, pos_x, pos_y, pos_z, densidades, N, Arbol, NodosParticulas)
		!write(*,*) "Arbol creado"
		detectados = 0
		lista_vecinos = 0
		distancias_vecinos = 0.0D+0
		distancias_vecinos_ord = 0.0D+0
		call Vecinos(NodosParticulas(idx), NodosParticulas, NodosParticulas(idx), N, detectados, n_vecinos, lista_vecinos, distancias_vecinos, 0, NIVEL_MAX_VECINOS, 0)
		!
		do while(.not. condicion_densidad)
			!Hay que acercar los vecinos a la particula, empezando de los más lejanos			
			do itr = 0, n_vecinos - 1, 1
				indices_vecinos(itr) = itr
			enddo

			distancias_vecinos_ord = ordenar(distancias_vecinos, indices_vecinos, n_vecinos)

			itr = n_vecinos - 1

			do while(itr >= 0 .and. .not. condicion_densidad)
				!write(*,*) pos_x(lista_vecinos(indices_vecinos(itr))), ",", pos_y(lista_vecinos(indices_vecinos(itr))), ",", pos_z(lista_vecinos(indices_vecinos(itr))), ",", distancias_vecinos_ord(itr)
				
				if(.not. visitados(lista_vecinos(indices_vecinos(itr)))) then
					visitados(lista_vecinos(indices_vecinos(itr))) = .true.
					write(*,*) "--__--", pos_x(lista_vecinos(indices_vecinos(itr))), ",", pos_y(lista_vecinos(indices_vecinos(itr))), ",", pos_z(lista_vecinos(indices_vecinos(itr))), ",0.,0.,0.,", distancias_vecinos_ord(itr), ",", masas(lista_vecinos(indices_vecinos(itr))), ",", densidad
				endif
				
				!Debo acercar la particula
				
				vector_posicion(0) = pos_x(lista_vecinos(indices_vecinos(itr)))
				vector_posicion(1) = pos_y(lista_vecinos(indices_vecinos(itr)))
				vector_posicion(2) = pos_z(lista_vecinos(indices_vecinos(itr)))
				
				!write(*,*) "ANTES: Posicion de vecino: ", vector_posicion, " Distancia: ", magnitudVector3D(diferenciaVectores3D(vector_posicion, vector_posicion_part))
				
				distancia = pos_x(lista_vecinos(indices_vecinos(itr))) - pos_x(idx)
				
				!write(*,*) "Distancia X: ", distancia, "---", pos_x(lista_vecinos(indices_vecinos(itr))), pos_x(idx)
				
				if(distancia > 0) then
					pos_x(lista_vecinos(indices_vecinos(itr))) = pos_x(lista_vecinos(indices_vecinos(itr))) - abs(distancia / 5.0D+0)
				else
					pos_x(lista_vecinos(indices_vecinos(itr))) = pos_x(lista_vecinos(indices_vecinos(itr))) + abs(distancia / 5.0D+0)
				endif

				distancia = pos_y(lista_vecinos(indices_vecinos(itr))) - pos_y(idx)

				!write(*,*) "Distancia Y: ", distancia, "---", pos_y(lista_vecinos(indices_vecinos(itr))), pos_y(idx)

				if(distancia > 0) then
					pos_y(lista_vecinos(indices_vecinos(itr))) = pos_y(lista_vecinos(indices_vecinos(itr))) - abs(distancia / 5.0D+0)
				else
					pos_y(lista_vecinos(indices_vecinos(itr))) = pos_y(lista_vecinos(indices_vecinos(itr))) + abs(distancia / 5.0D+0)
				endif

				distancia = pos_z(lista_vecinos(indices_vecinos(itr))) - pos_z(idx)

				!write(*,*) "Distancia Z: ", distancia, "---", pos_z(lista_vecinos(indices_vecinos(itr))), pos_z(idx)

				if(distancia > 0) then
					pos_z(lista_vecinos(indices_vecinos(itr))) = pos_z(lista_vecinos(indices_vecinos(itr))) - abs(distancia / 5.0D+0)
				else
					pos_z(lista_vecinos(indices_vecinos(itr))) = pos_z(lista_vecinos(indices_vecinos(itr))) + abs(distancia / 5.0D+0)
				endif
				
				vector_posicion(0) = pos_x(lista_vecinos(indices_vecinos(itr)))
				vector_posicion(1) = pos_y(lista_vecinos(indices_vecinos(itr)))
				vector_posicion(2) = pos_z(lista_vecinos(indices_vecinos(itr)))				
				
				!write(*,*) "DESPUES: Posicion de vecino: ", vector_posicion, " - Posicion particula: ", vector_posicion_part, " Distancia: ", magnitudVector3D(diferenciaVectores3D(vector_posicion, vector_posicion_part))
				
				!Compruebo si se alcanzó la densidad deseada
				
				densidad_local = calcularDensidadLocal(N, idx, lista_vecinos, n_vecinos, pos_x, pos_y, pos_z, masas)
				
				!write(*,*) "Densidad local para: ", idx, ": ", densidad_local
				
				if(densidad_local >= densidad * n_densidad) then
					write(*,*) "Logra condicion de densidad para particula: ", idx, densidad_local, " limite: ", densidad * n_densidad, " --- ", n_vecinos, " vecino: ", itr
					condicion_densidad = .true.
				endif
				
				itr = itr - 1
			enddo	
			
			!write(*,*) "Aun no logra condicion de densidad para: ", idx
		
		enddo
	enddo	
	
end subroutine 

!Calculo la velocidad inicial de cada particula basado en la velocidad angular w dada (radianes x seg) y la distancia de cada particula al eje
!de rotacion (z)
!
subroutine calcularDistribucionVelocidad(N, w, v_x, v_y, v_z, pos_x, pos_y, pos_z)

	use Vectores
	use Constantes
	
	implicit none
	
	integer N, i
	real*8 pos_x(0:N-1), pos_y(0:N-1), pos_z(0:N-1), v_x(0:N-1), v_y(0:N-1), v_z(0:N-1)
	real*8 vector_posicion2D(0:1), vector_posicion3D(0:2), vector_w(0:2), vector_velocidad(0:2)
	real*8 w, radio_xy, circunferencia, seccion

	do i = 0, N - 1, 1
		vector_posicion2D(0) = pos_x(i)
		vector_posicion2D(1) = pos_y(i)
		
		vector_posicion3D(0) = pos_x(i)
		vector_posicion3D(1) = pos_y(i)
		vector_posicion3D(2) = pos_z(i)
		
		radio_xy = magnitudVector2D(vector_posicion2D)
		circunferencia = 2.0D+0 * PI * radio_xy
		
		!Segun la velocidad angular en radianes, y sabiendo la "circunferencia" del punto cuyo centro es z
		!calculo la distancia recorrida por segundo (magnitud de la velocidad angular)
						
		seccion = (w * circunferencia) / (2.0D+0 * PI)
				
		!Ya puedo construir el vector de la velocidad angular (se rota sobre z)
		
		vector_w(0) = 0.0D+0
		vector_w(1) = 0.0D+0
		vector_w(2) = seccion
		
		!Con el producto cruz entre el vector posicion y el vector w, obtengo los componentes
		!de la velocidad de rotación
		
		vector_velocidad = productoCruz3D(vector_w, vector_posicion3D / magnitudVector3D(vector_posicion3D)) * PARSEC_MTS
		
		v_x(i) = vector_velocidad(0) 
		v_y(i) = vector_velocidad(1) 
		v_z(i) = vector_velocidad(2) 

!		write(*,*) "R ", radio_xy, "Circunferencia: ", circunferencia, "Distancia recorrida: ", seccion* PARSEC_MTS, " mts, Periodo: ", (2.0D+0 * PI) / w, " Velocidad rot: ", magnitudVector3D(vector_velocidad)
!		write(*,*) "----------"
		
		!write(*,*) i, "X,Y = ", pos_x(i), ",", pos_y(i), "Distancia X,Y de i = ", i, radio_xy, "Circun: ", circunferencia
	enddo

end subroutine

subroutine calcularPosicionesCilindro(N, Masa, Densidad, masas, pos_x, pos_y, pos_z, Radio_Nube, fraccion_altura)
	use Funciones
	use Constantes
	use Vectores
	use Fisica
	
	integer :: N, Ni, i, signo
	real*8 radio, Masa, Densidad, stop_cond, masa_densidad, randomsX, randomsY, randomsZ, Radio_Nube, fraccion_altura

	real*8 vector_posicion(0:1), vector_posicion_tmp(0:2)

	real*8 masas(0:N-1), pos_x(0:N-1), pos_y(0:N-1), pos_z(0:N-1)

	real*8 maxX, maxY, posZ

	CALL init_random_seed()
	
	radio = calcularRadioDisco(Masa, Densidad, fraccion_altura) / PARSEC_CMS		!en parsecs

	Radio_Nube = radio
	
	i = 0

	do while(i < N)
		stop_cond = Radio_Nube + 1.0D+0		

		!Esto es para asegurarse de que ningun punto exceda el radio local en que estoy trabajando
		do while(stop_cond > Radio_Nube)
			CALL RANDOM_NUMBER(randomsX)
			CALL RANDOM_NUMBER(randomsY)
			vector_posicion_tmp = asignarPosicion(Radio_Nube, randomsX, randomsY, -1.0D+0)
			vector_posicion(0) = vector_posicion_tmp(0)
			vector_posicion(1) = vector_posicion_tmp(1)
			stop_cond = magnitudVector2D(vector_posicion)
		enddo

		pos_x(i) = 0.0D+0 		
		pos_y(i) = vector_posicion(0)
		pos_z(i) = vector_posicion(1) 
		

		vector_posicion(0) = 0.0D+0
		vector_posicion(1) = 0.0D+0

		i = i + 1
	enddo

	do i = 0, N - 1, 1
		CALL RANDOM_NUMBER(randomsZ)
						
		posZ = frand(0.0D+0, Radio_Nube/fraccion_altura, randomsZ)
				
		pos_x(i) = posZ
	enddo

end subroutine

subroutine calcularPosicionesDisco(N, Masa, Densidad, masas, pos_x, pos_y, pos_z, Radio_Nube, fraccion_altura)
	use Funciones
	use Constantes
	use Vectores
	use Fisica
	
	integer :: N, Ni, i, signo
	real*8 radio, Masa, Densidad, stop_cond, masa_densidad, randomsX, randomsY, randomsZ, Radio_Nube, fraccion_altura

	real*8 vector_posicion(0:1), vector_posicion_tmp(0:2)

	real*8 masas(0:N-1), pos_x(0:N-1), pos_y(0:N-1), pos_z(0:N-1)

	real*8 maxX, maxY, posZ

	CALL init_random_seed()
	
	radio = calcularRadioDisco(Masa, Densidad, fraccion_altura) / PARSEC_CMS		!en parsecs

	Radio_Nube = radio
	
	i = 0

	do while(i < N)
		stop_cond = Radio_Nube + 1.0D+0		

		!Esto es para asegurarse de que ningun punto exceda el radio local en que estoy trabajando
		do while(stop_cond > Radio_Nube)
			CALL RANDOM_NUMBER(randomsX)
			CALL RANDOM_NUMBER(randomsY)
			vector_posicion_tmp = asignarPosicion(Radio_Nube, randomsX, randomsY, -1.0D+0)
			vector_posicion(0) = vector_posicion_tmp(0)
			vector_posicion(1) = vector_posicion_tmp(1)
			stop_cond = magnitudVector2D(vector_posicion)
		enddo

		pos_x(i) = vector_posicion(0)
		pos_y(i) = vector_posicion(1)
		pos_z(i) = 0.0D+0 

		vector_posicion(0) = 0.0D+0
		vector_posicion(1) = 0.0D+0

		i = i + 1
	enddo

	do i = 0, N - 1, 1
		CALL RANDOM_NUMBER(randomsZ)
						
		posZ = frand(0.0D+0, Radio_Nube/fraccion_altura, randomsZ)
				
		pos_z(i) = posZ
	enddo

end subroutine

subroutine calcularPosicionesEsfera(N, Masa, Densidad, masas, pos_x, pos_y, pos_z, Radio_Nube)
	use Constantes
	use Vectores
	use Fisica
	use Funciones
	
	integer :: N, Ni, i

	real*8 Masa, Densidad, stop_cond,  randomsX, randomsY, randomsZ, Radio_Nube

	real*8 vector_posicion(0:2)

	real*8 masas(0:N-1), pos_x(0:N-1), pos_y(0:N-1), pos_z(0:N-1)

	CALL init_random_seed()

	Radio_Nube = calcularRadio(Masa, Densidad) / PARSEC_CMS		!en parsecs

	i = 0

	do while(i < N)
		stop_cond = Radio_Nube + 1.0D+0		

		!Esto es para asegurarse de que ningun punto exceda el radio local en que estoy trabajando
		do while(stop_cond > Radio_Nube)
			CALL RANDOM_NUMBER(randomsX)
			CALL RANDOM_NUMBER(randomsY)
			CALL RANDOM_NUMBER(randomsZ)				
			vector_posicion = asignarPosicion(Radio_Nube, randomsX, randomsY, randomsZ)
			stop_cond = magnitudVector3D(vector_posicion)
		enddo

		pos_x(i) = vector_posicion(0)
		pos_y(i) = vector_posicion(1)
		pos_z(i) = vector_posicion(2)

		vector_posicion(0) = 0.0D+0
		vector_posicion(1) = 0.0D+0
		vector_posicion(2) = 0.0D+0
		i = i + 1
	enddo

end subroutine

!-----------------------------------------------
!---------     Programa Principal     ---------|
!-----------------------------------------------


program principal

	use Constantes
	use Funciones
	use Almacenamiento
	use Fisica
	
	implicit none

	!Masa en masas solares
	!Variacion en porcentaje de 0 a 1
	
	real*8 :: Masa_Nube, Masa_Nube_i, Variacion, Densidad_Nube, Radio_Nube, tff, w
	real*8 :: beta, altura, despl_x, despl_y, despl_z, n_densidad, veloc_x
	integer :: N, Ni, i, tipo, n_perturbacion	

	character(len=256) :: path
	
	character(len=100) :: cloud_title
	
	real*8 , allocatable, dimension(:) :: masas, pos_x, pos_y, pos_z, v_x, v_y, v_z, densidades
	character(256) :: namelistfile, prgname
	
	integer ipunit

	namelist /generateparam/ cloud_title,N,Masa_Nube,Densidad_Nube,Variacion,beta,tipo,despl_x,despl_y,despl_z,veloc_x,n_perturbacion,n_densidad,altura

	altura = 0.0D+0

	call getarg(0, prgname)
	call getarg(1, namelistfile)

	!write(*,*) "prgname y namelistfile - ", prgname, namelistfile

	ipunit = 100
	open(ipunit, file=namelistfile, status='old', action='read', err=100)
	read(ipunit, generateparam, err=104)
	close(ipunit)

	!print *, cloud_title,N,Masa_Nube,Densidad_Nube,Variacion,beta,tipo,despl_x,despl_y,despl_z,veloc_x,n_perturbacion,n_densidad,altura
	
	Masa_Nube_i = Masa_Nube
	Variacion = Variacion / 100.0

	Ni = N

	!----------------

	ALLOCATE(masas(0:N-1))
	ALLOCATE(pos_x(0:N-1))
	ALLOCATE(pos_y(0:N-1))
	ALLOCATE(pos_z(0:N-1))
	ALLOCATE(v_x(0:N-1))
	ALLOCATE(v_y(0:N-1))
	ALLOCATE(v_z(0:N-1))
	ALLOCATE(densidades(0:N-1))
	
	tff = freefalltimegcm(Densidad_Nube)/SEGS_YR

	do i = 0, N - 1, 1
		masas(i) = asignarMasa(Masa_Nube_i, Ni, Variacion)
		densidades(i) = Densidad_Nube
		!temporalmente, las velocidades empiezan en cero
		v_x(i) = 0.0D+0
		v_y(i) = 0.0D+0
		v_z(i) = 0.0D+0
	enddo
	
	if(tipo == 0) then
		call calcularPosicionesEsfera(N, Masa_Nube, Densidad_Nube, masas, pos_x, pos_y, pos_z, Radio_Nube)
	else if(tipo == 1) then
		call calcularPosicionesDisco(N, Masa_Nube, Densidad_Nube, masas, pos_x, pos_y, pos_z, Radio_Nube, altura)
	else if(tipo == 2) then
		call calcularPosicionesCilindro(N, Masa_Nube, Densidad_Nube, masas, pos_x, pos_y, pos_z, Radio_Nube, altura)
	else
		write(*,*) "Tipo de distribucion desconocida"
		return
	endif

	pos_x = pos_x + despl_x
	pos_y = pos_y + despl_y
	pos_z = pos_z + despl_z

	if(n_perturbacion > 0) then
		call PerturbarDensidad(N, pos_x, pos_y, pos_z, masas, Densidad_Nube, n_perturbacion, Radio_Nube, N/15, n_densidad)
	endif

	w = calculaVelocidadAngular(Masa_Nube*SOLAR_MASS_KG, Radio_Nube*PARSEC_MTS, beta)	
	
	if(.not. w == 0.0E+0) then
		call calcularDistribucionVelocidad(N, w, v_x, v_y, v_z, pos_x, pos_y, pos_z)
	endif
	
	write(*,*) 'Radio estimado: ', Radio_Nube, ' pc', " - ", Radio_Nube*PARSEC_CMS, " cms, Free fall time: ", tff, " Velocidad angular (rad/s): ", w, "- Masa Jeans para 10K: ", masaJeansH2(Densidad_Nube, 10.0D+0) / SOLAR_MASS_KG
	
	v_x = v_x + veloc_x
	
	path = "./datos/" // trim(cloud_title) // ".nc"

	call crearNubeNetCDF(path, N, Masa_Nube, Densidad_Nube, Variacion, beta, tipo, altura, Radio_Nube, tff, masaJeansH2(Densidad_Nube, 10.0D+0) / SOLAR_MASS_KG)
	call guardarNube(0, path, N, masas, pos_x, pos_y, pos_z, v_x, v_y, v_z, densidades)	
	call CerrarNubeCDF

	call actualizarSimParamsNetCDF(path, 10.0D+0, 50.0D+0, 0.7D+0, 35, 100)
	call CerrarNubeCDF
	stop

100	write (6, * ) 'Cannot read namelist: generateparam (1)', trim(namelistfile)
	stop
104	write ( 6, * ) 'Cannot read namelist: generateparam (2)',  trim(namelistfile)

end
