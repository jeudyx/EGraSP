MODULE BarnesHut
	CONTAINS

!Calcula la aceleración gravitacional que el conjunto de particulas dentro del arbol (octree) ejerce sobre la particula
!siguiendo el algoritmo Barnes-Hut http://arborjs.org/docs/barnes-hut

!Si umbral es cero, entonces se convierte en fuerza bruta

!Resultado en m/s^2 ?

recursive function calcularAccGravBH(p, Arbol, umbral, soft_len) result(acc_vect)

	use Constantes
	use Auxiliar
	use Tipos
	use Vectores
	use Funciones
	use Octree
	use Fisica

	implicit none

	real*8 acc_vect(0:2), tmp_acc_vect(0:2), centro_masa(0:2)
	real*8 umbral, lado, distancia12, soft_len, radio1, radio2
	type(Particula) p, q

	type(OctreeNode), POINTER :: Arbol, Hijo

	type(Cubo) cube

	logical son_hermanos

	integer i

	acc_vect(0) = 0.0D+0
	acc_vect(1) = 0.0D+0
	acc_vect(2) = 0.0D+0	

	!Si no es hoja, y no tiene particulas, es un nodo vacio, no hacer nada
	if(.not. Arbol%hoja .and. Arbol%n_particulas == 0) then
		return
	endif

	cube = CrearCubo(Arbol%radio, Arbol%centroide)	

	if(Arbol%hoja) then	
		!que no sea la misma particula
		if(.not. Arbol%id_particula == p%id .and. Arbol%id_particula > 0) then
!			write(*,*) "Llego a hoja. Particula: ", Arbol%id_particula
			q = construirParticula(Arbol%id_particula, Arbol%masa/SOLAR_MASS_KG, Arbol%centro_masa(0), Arbol%centro_masa(1), Arbol%centro_masa(2), 0.0D+0, 0.0D+0, 0.0D+0, Arbol%densidad)

			!write(*,*) "Creando softening leght. Densidades: ", p%densidad, q%densidad, q%masa, q%id

			acc_vect = aceleracion_g_vect(p%posicion*PARSEC_MTS, q%posicion*PARSEC_MTS, q%masa, soft_len)
			return
		else
			!devuelve cero
			return
		endif
	else
		!El lado del cubo es la distancia entre el vertice 1 y el 2 (distancia en metros)

		lado = distanciaPuntos(cube%vertice1, cube%vertice2) * PARSEC_MTS

		!La comparacion con el umbral es lado / distancia

		q = construirParticula(Arbol%id_particula, Arbol%masa/SOLAR_MASS_KG, Arbol%centro_masa(0), Arbol%centro_masa(1), Arbol%centro_masa(2), 0.0D+0, 0.0D+0, 0.0D+0, 0.0D+0)

		!write(*,*) "Comparando posiciones: P: ", p%posicion, "Q: ", q%posicion, "Distancia cruda: ", distanciaParticulas(p,q)

		distancia12 = distanciaParticulas(p,q)


		if( (lado / distancia12) < umbral) then

		 	!If s/d < umbral, treat this internal node as a single body, and calculate the force it exerts on body b, and add this amount to b’s net force.

		 	acc_vect = aceleracion_g_vect(p%posicion*PARSEC_MTS, q%posicion*PARSEC_MTS, q%masa, soft_len)

		 	return
		else
			!Otherwise, run the procedure recursively on each of the current node’s children

!			write(*,*) "s/d >= umbral. Examinando hijos."

			do i = 0, 7, 1
				tmp_acc_vect(0) = 0.0D+0
				tmp_acc_vect(1) = 0.0D+0
				tmp_acc_vect(2) = 0.0D+0

!				write(*,*) "Sobre Hijo: ", i, " ID nodo: ", Arbol%hijos(i)%id

				Hijo => Arbol%hijos(i)

				tmp_acc_vect = calcularAccGravBH(p, Hijo, umbral, soft_len)
				acc_vect(0) = acc_vect(0) + tmp_acc_vect(0)
				acc_vect(1) = acc_vect(1) + tmp_acc_vect(1)
				acc_vect(2) = acc_vect(2) + tmp_acc_vect(2)
			enddo
!			write(*,*) "------"
			return
		endif
	endif
	
	
	return

end function

END MODULE BarnesHut
