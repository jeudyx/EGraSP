MODULE Dinamica
	CONTAINS

!Procedimiento para avanzar el estado dinámico del sistema un timestep, utilizando el método Leap-Frog

!Este será llamado desde afuera por un numero de iteraciones deseado. 
!En la primera iteración, se debe calcular los v a partir de la aceleracion inicialmente calculada
!Recordar que aceleraciones se definen a mitad del intervalo de tiempo dt, la posición y aceleracion en el punto mismo, por
!eso en la primera iteracion se calcula la mitad, y la otra mitad se calcula en la presente iteracion, etc

!Se asume que en acc_i vienen los componentes de las aceleraciones para todas las particulas, previamente calculadas
!con Barnes-Hut y otros componentes (presion, campo magnético, etc).

!El dt viene en años, hay que traducir a segundos

!Al 10/12/2011 se está usando el esquema leapfrog descrito en el paper 1 de Gadget, paginas 12 en adelante

subroutine pasoLeapFrog(N, itr_inicio, itr_final, Arbol, NodosParticulas, masas, pos_x, pos_y, pos_z, densidades, densidades_locales, v_x, v_y, v_z, acc_x, acc_y, acc_z, dt, umbralBH, tolerancia_colision, beta, n_vecinos, matriz_vecinos, presiones, soft_len, temperatura, myid)

	use Constantes
	use Tipos
	use Funciones
	use Fuerzas
	use Octree
	use BarnesHut
	use Vectores
	use Fisica
	
	implicit none

	integer N, i, j, itr_inicio, itr_final, myid, n_vecinos, max_i
	integer lista_vecinos(0:n_vecinos-1)
	integer matriz_vecinos(0:N-1, 0:n_vecinos-1)
	real*8 presiones(0:N-1)
	real*8 pos_x(0:N-1), pos_y(0:N-1), pos_z(0:N-1), acc_x(0:N-1), acc_y(0:N-1), acc_z(0:N-1), masas(0:N-1), densidades(0:N-1), densidades_locales(0:N-1), pos_predictor_x(0:N-1), pos_predictor_y(0:N-1), pos_predictor_z(0:N-1)
	real*8 v_x(0:N-1), v_y(0:N-1), v_z(0:N-1), velocidades_angulares(0:N-1)
	real*8 dt, umbralBH, dist_max, dist_avg, tolerancia_colision, beta, temperatura, densidad, soft_len, mag
	real*8 acc_grav_vect(0:2), vector_posicion(0:2), vector_velocidad(0:2), acc_presion_vect(0:2), acc_visc_vect(0:2), gradiente_presion(0:2), fuerza_presion_vect(0:2)

	type(Particula) p

	type(OctreeNode), POINTER :: Arbol
	
	type(OctreeNode) NodosParticulas(0:N-1)
	
	dist_max = 0.0D+0
	dist_avg = 0.0D+0

!	write(*,*) myid, "pasoLeapFrog, beta: ", beta

!	write(*,*) myid, "Entrando en pasoLeapFrog. Posiciones(x): ", pos_x, " - \nPosiciones(y): ", pos_y, "- \nPosiciones(z): ", pos_z

	!Temporalmente, predigo la posicion en el siguiente medio timestep, para actualizar la aceleraci�n
	do i = 0, N - 1, 1
		!
		!En este punto, las velocidades están calculadas para el step k-1/2 
		!(k siendo el timestep actual desde el ciclo externo a esta funcion desde donde se llama		
		pos_predictor_x(i) = ((pos_x(i) * PARSEC_MTS) + (v_x(i) * (dt * SEGS_YR * 0.5D+0))) / PARSEC_MTS
		pos_predictor_y(i) = ((pos_y(i) * PARSEC_MTS) + (v_y(i) * (dt * SEGS_YR * 0.5D+0))) / PARSEC_MTS
		pos_predictor_z(i) = ((pos_z(i) * PARSEC_MTS) + (v_z(i) * (dt * SEGS_YR * 0.5D+0))) / PARSEC_MTS
		!	
	
		!Necesito saber el radio/distancia maxima para reconstruir arbol
		vector_posicion(0) = pos_predictor_x(i)
		vector_posicion(1) = pos_predictor_y(i)
		vector_posicion(2) = pos_predictor_z(i)
		
		mag = magnitudVector3D(vector_posicion)

		dist_avg = dist_avg + mag
		
		if(dist_max < mag) then			
			dist_max = mag
			max_i = i
		endif		
		
	enddo	
	
	dist_avg = dist_avg / N
	
	!En cada caso se debe reconstruir el arbol
	!uso las nuevas posiciones. Esta correccion se hace el 10/08 por error de particula escapando

	call limpiarArbol(Arbol)
	Arbol%radio = dist_max	
	!write(*,*) myid, "-Maxima distancia: ", dist_max, " en particula: ", max_i, "posiciones: ", pos_predictor_x(max_i), pos_predictor_y(max_i), pos_predictor_z(max_i), "Velocidad: ", v_x(max_i), v_y(max_i), v_z(max_i)
	call CrearOctree(masas, pos_predictor_x, pos_predictor_y, pos_predictor_z, densidades, N, Arbol, NodosParticulas)

	!write(*,*) "Octree reconstruido en leapfrog step"

	do i = itr_inicio, itr_final, 1
		
		vector_posicion(0) = pos_predictor_x(i)
		vector_posicion(1) = pos_predictor_y(i)
		vector_posicion(2) = pos_predictor_z(i)
		
		if(magnitudVector3D(vector_posicion) >= FACTOR_DISTANCIA_MAX * dist_avg) then
			!write(*,*) i, " - ALERTA: posicion de particula: ", vector_posicion, " esta a una distancia ", magnitudVector3D(vector_posicion)," superior al limite: ", FACTOR_DISTANCIA_MAX * dist_avg, ". Distancia promedio: ", dist_avg
			v_x(i) = 0.0D+0
			v_y(i) = 0.0D+0
			v_z(i) = 0.0D+0
		else
		
			!Se calcula la velocidad del siguiente step. Para esto se necesita recalcular la aceleración
			!con esto queda la aceleración y velocidad lista para el siguiente timestep

			!Utilizo las posiciones y velocidades actuales para el calculo de la aceleracion

			acc_grav_vect(:) = 0.0D+0
			acc_presion_vect(:) = 0.0D+0
			acc_visc_vect(:) =  0.0D+0


			!Calculo la aceleración basado en las posiciones predichas
			p = construirParticula(i, masas(i), pos_predictor_x(i), pos_predictor_y(i), pos_predictor_z(i), v_x(i), v_y(i), v_z(i), densidades(i))
			acc_grav_vect = calcularAceleracion(p, Arbol, umbralBH, soft_len, N)

			!Estimo el gradiente de presión, y la aceleracion

			lista_vecinos = matriz_vecinos(i, 0:n_vecinos-1)

			gradiente_presion = GradientePresion(i, n_vecinos, lista_vecinos, presiones, N, pos_predictor_x, pos_predictor_y, pos_predictor_z, masas, densidades, soft_len)

			acc_presion_vect = gradiente_presion / (densidades_locales(i) / 1000.0D+0)

			!Ecuacion de movimiento: gravedad - presion
				
			acc_x(i) = acc_grav_vect(0) - acc_presion_vect(0)
			acc_y(i) = acc_grav_vect(1) - acc_presion_vect(1)
			acc_z(i) = acc_grav_vect(2) - acc_presion_vect(2)
	!!		write(*,*) "-------------------------------------"

			vector_posicion(0) = pos_predictor_x(i)
			vector_posicion(1) = pos_predictor_y(i)
			vector_posicion(2) = pos_predictor_z(i)

			if(beta > 0.0D+0) then			
				velocidades_angulares(i) = velocidad_angular(acc_grav_vect, vector_posicion * PARSEC_MTS, beta)
				vector_velocidad = vector_velocidad_angular(vector_posicion, velocidades_angulares(i))
			else
				velocidades_angulares(i) = 0.0D+0
				vector_velocidad = 0.0D+0
			endif
			
			v_x(i) = v_x(i) + (acc_x(i) * (dt * SEGS_YR))
			v_y(i) = v_y(i) + (acc_y(i) * (dt * SEGS_YR))
			v_z(i) = v_z(i) + (acc_z(i) * (dt * SEGS_YR))

			!01/16/2012 Cambio manera de incorporar viscosidad artificial
			!Con la velocidad predicha, calculo la viscosidad artificial

!!!			acc_visc_vect = ArtificialViscosityAcc(i, temperatura, soft_len, n_vecinos, lista_vecinos, N, pos_predictor_x, pos_predictor_y, pos_predictor_z, v_x, v_y, v_z, densidades, densidades_locales, masas)

			if(.false. .and. (p%id == 100 .or. p%id == 500 .or. p%id == 990)) then
				!write(*,*) myid, i, " En dinamica.pasoLeapFrog. Particula: ", p%id, " - aceleracion neta sin visc: ", acc_x(i), acc_y(i), acc_z(i), " - aceleracion gravedad: ", acc_grav_vect, magnitudVector3D(acc_grav_vect), " - aceleracion presion: ", acc_presion_vect, magnitudVector3D(acc_presion_vect), " Aceleracion viscosidad: ", acc_visc_vect, magnitudVector3D(acc_visc_vect), " - Gradiente presion: ", magnitudVector3D(gradiente_presion), " - Velocidades: ", v_x(i), v_y(i), v_z(i), " - Densidad(i): ", densidades_locales(i)
				write(*,*) p%id, " - aceleracion gravedad: ", acc_grav_vect, magnitudVector3D(acc_grav_vect), " - aceleracion presion: ", acc_presion_vect, magnitudVector3D(acc_presion_vect), "- aceleracion neta sin visc: ", acc_x(i), acc_y(i), acc_z(i), " Aceleracion viscosidad: ", acc_visc_vect, magnitudVector3D(acc_visc_vect), " - Gradiente presion: ", gradiente_presion, magnitudVector3D(gradiente_presion), " - Densidad(i): ", densidades_locales(i)
			endif

			!Incorporo componentes de velocidad por aceleracion de viscosidad artificial

!!!			v_x(i) = v_x(i) + (acc_visc_vect(0) * (dt * SEGS_YR))
!!!			v_y(i) = v_y(i) + (acc_visc_vect(1) * (dt * SEGS_YR))
!!!			v_z(i) = v_z(i) + (acc_visc_vect(2) * (dt * SEGS_YR))

			!Actualizo las posiciones reales con esta velocidad (del siguiente medio ts) y la actual

			pos_x(i) = ((pos_x(i) * PARSEC_MTS) + ( (v_x(i) + vector_velocidad(0)) * (dt * SEGS_YR * 0.5D+0))) / PARSEC_MTS
			pos_y(i) = ((pos_y(i) * PARSEC_MTS) + ( (v_y(i) + vector_velocidad(1)) * (dt * SEGS_YR * 0.5D+0))) / PARSEC_MTS
			pos_z(i) = ((pos_z(i) * PARSEC_MTS) + ( (v_z(i) + vector_velocidad(2)) * (dt * SEGS_YR * 0.5D+0))) / PARSEC_MTS
		endif			
	enddo

	!write(*,*) "Aceleraciones actualizadas"

end subroutine

SUBROUTINE calcularVecinos_Densidad_Presion(N, itr_inicio, itr_final, Arbol, NodosParticulas, masas, pos_x, pos_y, pos_z, v_x, v_y, v_z, n_vecinos, densidades_locales, matriz_vecinos, presiones, temperatura, myid)

	use Constantes
	use Tipos
	use Octree
	use Vectores
	use Fisica
	
	implicit none

	real*8 temperatura

	real*8 pos_x(0:N-1), pos_y(0:N-1), pos_z(0:N-1), masas(0:N-1), v_x(0:N-1), v_y(0:N-1), v_z(0:N-1), densidades_locales(0:N-1), presiones(0:N-1), distancias_vecinos(0:n_vecinos-1)
	
	integer N, i, itr_inicio, itr_final, n_vecinos, detectados, myid, j

	integer lista_vecinos(0:n_vecinos - 1)

	integer matriz_vecinos(0:N-1,0:n_vecinos - 1)

	type(Particula) p

	type(OctreeNode), POINTER :: Arbol
	
	type(OctreeNode) NodosParticulas(0:N-1)
	
	do i = itr_inicio, itr_final, 1	
		detectados = 0
		lista_vecinos = 0
		distancias_vecinos = 0.0D+0
		
		call ResetearVisitados(Arbol)
		
		do j = 0, N - 1, 1
			NodosParticulas(j)%visitado = .false.
		enddo
				
		call Vecinos(NodosParticulas(i), NodosParticulas, NodosParticulas(i), N, detectados, n_vecinos, lista_vecinos, distancias_vecinos, 0, NIVEL_MAX_VECINOS, myid)
		matriz_vecinos(i, 0:n_vecinos-1) = lista_vecinos

		densidades_locales(i) = calcularDensidadLocal(N, i, lista_vecinos, n_vecinos, pos_x, pos_y, pos_z, masas)
		presiones(i) = PresionGasIdeal(densidades_locales(i), temperatura)

		if(lista_vecinos(n_vecinos - 1) == 0 .and. lista_vecinos(n_vecinos - 2) == 0 .and. lista_vecinos(n_vecinos - 3) == 0) then
			write(*,*) "ALERTA: vecinos nullos encontrados para particula: ", i, " Densidad local: ", densidades_locales(i), " - Presion: ", presiones(i), " - Vecinos: ", lista_vecinos
		endif
				
	enddo
	
END SUBROUTINE

END MODULE Dinamica

