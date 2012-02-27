MODULE BarnesHut
	CONTAINS

!Calcula la aceleración gravitacional que el conjunto de particulas dentro del arbol (octree) ejerce sobre la particula
!siguiendo el algoritmo Barnes-Hut http://arborjs.org/docs/barnes-hut

!Si umbral es cero, entonces se convierte en fuerza bruta

!Resultado en m/s^2 ?

recursive function calcularAccGravBH(p, Arbol, umbral) result(acc_vect)

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

	!write(*,*) "Umbral BH: ", umbral

!	write(*,*) "------"
!	write(*,*) "Entro en calcularAccGravBH. Arbol (id): ", Arbol%id, " hoja? ", Arbol%hoja, " Particulas: ", Arbol%n_particulas, "Centro masa: ", Arbol%centro_masa, " Masa: ", Arbol%masa / SOLAR_MASS_KG
	!write(*,*) "Posicion de la particula: ", p%posicion, " Masa: ", p%masa / SOLAR_MASS_KG

	cube = CrearCubo(Arbol%radio, Arbol%centroide)	

	if(Arbol%hoja) then	
		!que no sea la misma particula
		if(.not. Arbol%id_particula == p%id .and. Arbol%id_particula > 0) then
!			write(*,*) "Llego a hoja. Particula: ", Arbol%id_particula
			q = construirParticula(Arbol%id_particula, Arbol%masa/SOLAR_MASS_KG, Arbol%centro_masa(0), Arbol%centro_masa(1), Arbol%centro_masa(2), 0.0D+0, 0.0D+0, 0.0D+0, Arbol%densidad)
			
			!write(*,*) "Creando softening leght. Densidades: ", p%densidad, q%densidad, q%masa, q%id
			
			radio1 = calcularRadio(p%masa / SOLAR_MASS_KG, p%densidad) / PARSEC_CMS
			radio2 = calcularRadio(q%masa / SOLAR_MASS_KG, q%densidad) / PARSEC_CMS
						
			soft_len = (radio1 + radio2) * PARSEC_MTS				
			!soft_len = soft_len * 1.5
													
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

		radio1 = calcularRadio(p%masa / SOLAR_MASS_KG, p%densidad) / PARSEC_CMS

		q = construirParticula(Arbol%id_particula, Arbol%masa/SOLAR_MASS_KG, Arbol%centro_masa(0), Arbol%centro_masa(1), Arbol%centro_masa(2), 0.0D+0, 0.0D+0, 0.0D+0, 0.0D+0)
		
		!write(*,*) "Comparando posiciones: P: ", p%posicion, "Q: ", q%posicion, "Distancia cruda: ", distanciaParticulas(p,q)
		
		distancia12 = distanciaParticulas(p,q)
		
!		write(*,*) "Lado: ", lado / PARSEC_MTS, " distancia: ", distancia12/ PARSEC_MTS, " s/d: ", lado / distancia12
		
		if( (lado / distancia12) < umbral) then
		 	
!		 	write(*,*) "s/d < umbral. Este nodo interno se tratará como particula: Nodo: ", Arbol%id
		 	
		 	!If s/d < umbral, treat this internal node as a single body, and calculate the force it exerts on body b, and add this amount to b’s net force.
		 	
			soft_len = (radio1 + radio1) * PARSEC_MTS			
			!soft_len = soft_len * 1.5		 	
			
		 	
		 	acc_vect = aceleracion_g_vect(p%posicion*PARSEC_MTS, q%posicion*PARSEC_MTS, q%masa, soft_len)
			
!			if(magnitudVector3D(acc_vect) >= 1.0E-7) then
!				write(*,*) p%id, p%masa/SOLAR_MASS_KG, "- Posicion particula 57: ", p%posicion/PARSEC_MTS, " - Aceleracion BH nodo interno para particula : ", acc_vect, " Datos nodo: posicion: ", q%posicion/PARSEC_MTS, " - masa: ", q%masa/SOLAR_MASS_KG, " Distancia: ", distancia12/ PARSEC_MTS, " - Radio particula: ", radio1, " lado: ", lado
!			endif
			
!		 	write(*,*) "------ Resultado de aceleracion con nodo interno cercano: ", acc_vect
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
				tmp_acc_vect = calcularAccGravBH(p, Hijo, umbral)
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
