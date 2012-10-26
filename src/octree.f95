MODULE Octree
	CONTAINS

RECURSIVE SUBROUTINE imprimirArbol(Arbol)

	use Tipos
	
	implicit none

	integer*4 i
		
	type(OctreeNode), POINTER :: Arbol, Hijo

	write(*,*) "Nodo: ", Arbol%id
	write(*,*) "----"

	if(.not. Arbol%hoja .and. Arbol%hijos_creados) then
		do i = 0, 7, 1
			Hijo => DarHijo(Arbol, i)
			write(*,*) "Hijo: ", Hijo%id, Hijo%hoja, Hijo%hijos_creados
		enddo
		write(*,*) "--------"
		write(*,*) "--------"	
		do i = 0, 7, 1
			Hijo => DarHijo(Arbol, i)
			call imprimirArbol(Hijo)
		enddo
	endif
	
END SUBROUTINE


RECURSIVE SUBROUTINE limpiarArbol(Arbol)

	use Tipos
	
	implicit none

	integer*4 i
		
	type(OctreeNode), POINTER :: Arbol, Hijo

!	write(*,*) "Limpiado nodo: ", Arbol%id

	if(.not. Arbol%hoja .and. Arbol%hijos_creados) then

		do i = 0, 7, 1
			Hijo => DarHijo(Arbol, i)
			call limpiarArbol(Hijo)
		enddo
	endif
	
!	write(*,*) "Deallocating nodo: ", Arbol%id
	
!	write(*,*) "--------"
!	write(*,*) "--------"	

	deallocate(Arbol)

END SUBROUTINE


RECURSIVE SUBROUTINE xlimpiarArbol(Arbol)

	use Tipos
	
	implicit none

	integer*4 i
		
	type(OctreeNode), POINTER :: Arbol, Hijo
	
	Arbol%id = 0
	Arbol%id_particula = -1
	Arbol%n_particulas = 0
	Arbol%hoja = .true.
	Arbol%centroide(0) = 0.0D+0
	Arbol%centroide(1) = 0.0D+0
	Arbol%centroide(2) = 0.0D+0
	Arbol%centro_masa(0) = 0.0D+0
	Arbol%centro_masa(1) = 0.0D+0
	Arbol%centro_masa(2) = 0.0D+0
	Arbol%radio = 0.0D+0
	Arbol%masa = 0.0D+0
	Arbol%densidad = 0.0D+0			
		
	if(Arbol%hijos_creados) then	
		do i = 0, 7, 1
			Hijo => DarHijo(Arbol, i)
			Hijo%id = 0
			Hijo%id_particula = -1
			Hijo%n_particulas = 0
			Hijo%hoja = .true.
			Hijo%centroide(0) = 0.0D+0
			Hijo%centroide(1) = 0.0D+0
			Hijo%centroide(2) = 0.0D+0
			Hijo%centro_masa(0) = 0.0D+0
			Hijo%centro_masa(1) = 0.0D+0
			Hijo%centro_masa(2) = 0.0D+0
			Hijo%radio = 0.0D+0
			Hijo%masa = 0.0D+0
			Hijo%densidad = 0.0D+0
			call limpiarArbol(Hijo)
		enddo
	endif

END SUBROUTINE

!Genera un octree para calcular fuerzas mediante el algoritmo Barnes-Hut descrito en http://arborjs.org/docs/barnes-hut

SUBROUTINE CrearOctree(masas, coordenadas_x, coordenadas_y, coordenadas_z, densidades, N, Arbol, NodosParticulas)

	use Constantes
	use Tipos
	use Funciones
	
	implicit none

	integer*4 N, i, NNodos, k
	
	real*8 masas(0:N-1), coordenadas_x(0:N-1), coordenadas_y(0:N-1), coordenadas_z(0:N-1), densidades(0:N-1)

	type(OctreeNode) NodosParticulas(0:N-1)

	type(OctreeNode), POINTER :: Arbol	!Cabeza del arbol
	type(Particula) p
	
	type(Cubo) cube

	real*8 posicion(0:2)
		
	!Teniendo el árbol en memoria, inserto partícula a partícula según el algoritmo descrito en http://arborjs.org/docs/barnes-hut
		
	!Cubo que contiene la nube entera	
	
	!Centro de la nube
	posicion(0) = 0.0D+0
	posicion(1) = 0.0D+0
	posicion(2) = 0.0D+0

	cube = CrearCubo(Arbol%radio, posicion)
	
	NNodos = 1
	
	Arbol%nivel = 0
	
	!Inserta las partículas una a una
	do i = 0, N - 1, 1
	
		!Es un esquema top-down, cada particula la intentará insertar siempre desde la raíz		
		
		if(masas(i) > 0.0D+0) then
			!Solo inserto particulas existentes, no eliminadas
			p = construirParticula(i, masas(i), coordenadas_x(i), coordenadas_y(i), coordenadas_z(i), 0.0D+0, 0.0D+0, 0.0D+0, densidades(i))
			!write(*,*) "Insertando particula: ", i, " a arbol. NNodos = ", NNodos
			call InsertarParticula(p, Arbol, NNodos, masas, coordenadas_x, coordenadas_y, coordenadas_z, N, NodosParticulas)
		endif
	enddo

END SUBROUTINE

!Inserta una par	ticula en el arbol. "i" es la posicion del nodo dentro del árbol donde se intentará meter la partícula
!k es un apuntador de los nodos usados hasta el momento. Apunta al siguiente nodo disponible para ser usado en el arbol

RECURSIVE SUBROUTINE InsertarParticula(p, Arbol, NNodos, masas, coordenadas_x, coordenadas_y, coordenadas_z, N, NodosParticulas)
	use Tipos
	use Funciones
	use Fisica
	implicit none

	integer*4 NNodos, j,  N, id_padre
	real*8 masas(0:N-1), coordenadas_x(0:N-1), coordenadas_y(0:N-1), coordenadas_z(0:N-1)
	real*8 tmp_centromasas(0:2)
	

	type(OctreeNode), POINTER ::  Hijo, ArbolPtr

	type(OctreeNode), TARGET ::  Arbol
	type(Particula) p, q
	type(Cubo) hijos_nuevos(0:7)
	type(Cubo) cube, tmp_cube
	
	type(OctreeNode) NodosParticulas(0:N-1)
	
	logical encontrado_nuevo, encontrado_existente
	
	encontrado_nuevo = .false.
	encontrado_existente = .false.
	
	ArbolPtr => Arbol
	
	cube = CrearCubo(Arbol%radio, Arbol%centroide)
	
	if(Arbol%id_particula == -1) then
	
		!Podria ser un nodo interno o bien vacio
	
		if(Arbol%n_particulas == 0) then	
			!write(*,*) "Nodo ", Arbol%id, " es un nodo vacio, asigno particula aqui: ", p%id, " hoja? ", Arbol%hoja
			!Nodo vacio
			!If node x does not contain a body, put the new body p here.
			Arbol%id_particula = p%id
			Arbol%n_particulas = 1
			Arbol%masa = p%masa
			Arbol%densidad = p%densidad
			Arbol%centro_masa(0) = p%posicion(0)
			Arbol%centro_masa(1) = p%posicion(1)
			Arbol%centro_masa(2) = p%posicion(2)
			Arbol%hoja = .true. 
			
			NodosParticulas(p%id) = Arbol
						
		else
			!Nodo interno (contiene particulas en subnodos)
			!OJO, asume que los hijos estan inicializados
						
			!If node is an internal node, update the center-of-mass and total mass of node. 
			!Recursively insert the body b in the appropriate quadrant.

			Arbol%n_particulas = Arbol%n_particulas + 1
						
			tmp_centromasas = centroMasa(p%masa, p%posicion(0), p%posicion(1), p%posicion(2), Arbol%masa, Arbol%centro_masa(0), Arbol%centro_masa(1), Arbol%centro_masa(2))
			
			Arbol%masa = Arbol%masa + p%masa
			
			Arbol%centro_masa(0) = tmp_centromasas(0)
			Arbol%centro_masa(1) = tmp_centromasas(1)
			Arbol%centro_masa(2) = tmp_centromasas(2)
			
			Arbol%hoja = .false.
			
			do j = 0, 7, 1
				!Busco cual hijo contiene la particula
				!obtengo el cubo que corresponde a cada hijo
				Hijo => DarHijo(ArbolPtr, j)
				tmp_cube = CrearCubo(Hijo%radio, Hijo%centroide)
				!Si el cubo hijo contiene a la particula existente
				if(CuboContienePunto(tmp_cube, p%posicion)) then				
					call InsertarParticula(p, Hijo, NNodos, masas, coordenadas_x, coordenadas_y, coordenadas_z, N, NodosParticulas)
					return
				endif				
			enddo						
			write(*,*) "OJO, en nodo interno, nunca encontró hijo. Particula: ", p%id, "Nodo: (id, id particula, es hoja?, n_particulas) ", Arbol%id, Arbol%id_particula, Arbol%hoja, Arbol%n_particulas,  "Posicion particula: ", p%posicion
						
			!Debugging de problema de no deteccion en 10/01/2011
			if(.true.) then

				write(*,*) "--------------------------"

				call imprimirNube(N, masas, coordenadas_x, coordenadas_y, coordenadas_z, coordenadas_x, coordenadas_y, coordenadas_z)

				write(*,*) "--------------------------"

				write(*,*) "Cubo de nodo actual"
				tmp_cube = CrearCubo(Arbol%radio, Arbol%centroide)
				call imprimirCubo(tmp_cube)
				write(*,*) "Cubo dado como argumento"
				call imprimirCubo(cube)
				write(*,*) "---------------------------------"

				j = 0
				!mato al programa para analizar
				j = 10 / j
				write(*,*) "MUERE"				
				call MPI_abort (1)
			endif
		endif
	else			
		!Nodo externo: era una hoja, se tendrá que dividir
					
		!If node x is an external node, say containing a body named c, then there are two bodies b and c in the same region
		!Subdivide the region further by creating four children. 
		
		!Divido el cubo en 8 subregiones, cada una se le asignará a un hijo, tomando para cada hijo
		!un nodo del "pool" (arbol), avanzando k


		!Es el root, no tiene padre
		if(Arbol%id == 0) then
			id_padre = -1
		else
			id_padre = Arbol%padre%id
		endif
		
		!write(*,*) "Metiendo a nodo externo lleno: (i, id_part existente, id part new)", Arbol%id, Arbol%id_particula, p%id, " Cubo a dividir. Tiene hijos: ", Arbol%hijos_creados, " Padre: ", id_padre
		
		
		!Creo cubos a partir del cubo padre
		hijos_nuevos = DividirCubo(cube)
		!Inicializo cada nodo
		
		!Creo espacio para los hijos si no existia ya
		if(.not. Arbol%hijos_creados) then
			allocate (Arbol%hijo0)
			allocate (Arbol%hijo1)
			allocate (Arbol%hijo2)
			allocate (Arbol%hijo3)
			allocate (Arbol%hijo4)
			allocate (Arbol%hijo5)
			allocate (Arbol%hijo6)
			allocate (Arbol%hijo7)

			Arbol%hijo0%hijos_creados = .false.
			Arbol%hijo1%hijos_creados = .false.
			Arbol%hijo2%hijos_creados = .false.
			Arbol%hijo3%hijos_creados = .false.
			Arbol%hijo4%hijos_creados = .false.
			Arbol%hijo5%hijos_creados = .false.
			Arbol%hijo6%hijos_creados = .false.
			Arbol%hijo7%hijos_creados = .false.

			Arbol%hijos_creados = .true.
		endif
						
		!Inicializo/Reseteo los hijos
		
		Arbol%hoja = .false.
				
		do j = 0, 7, 1
			Hijo => DarHijo(ArbolPtr, j)
			Hijo%id = NNodos
			Hijo%id_particula = -1
			Hijo%hoja = .true. 
			if(.not. Arbol%hijos_creados) then
				!Si es la primera vez que el hijo se crea, se asigna, sino, se mantiene el valor anterior
				Hijo%hijos_creados = .false. 
			endif
			
			Hijo%n_particulas = 0
			Hijo%centroide(0) = hijos_nuevos(j)%centroide(0)
			Hijo%centroide(1) = hijos_nuevos(j)%centroide(1)
			Hijo%centroide(2) = hijos_nuevos(j)%centroide(2)
			Hijo%centro_masa(0) = -1.0D+0
			Hijo%centro_masa(1) = -1.0D+0
			Hijo%centro_masa(2) = -1.0D+0
			Hijo%radio = hijos_nuevos(j)%radio
			Hijo%padre => Arbol
			Hijo%nivel = Arbol%nivel + 1
			NNodos = NNodos + 1
		enddo				

		if(.not. Arbol%hijos_creados) then
			
			!Voy a asignar los vecinos a los hijos recien creados
			do j = 0, 7, 1
				call AsignarAdyacentes(j, Hijo, Arbol)
			enddo
			
			Arbol%hijos_creados = .true.
		endif

					
		!Then, recursively insert both b and c into the appropriate quadrant(s). 
		!Since b and c may still end up in the same quadrant,
		!there may be several subdivisions during a single insertion. Finally, update the center-of-mass and total mass of x


		Arbol%n_particulas = Arbol%n_particulas + 1

		tmp_centromasas = centroMasa(p%masa, p%posicion(0), p%posicion(1), p%posicion(2), Arbol%masa, Arbol%centro_masa(0), Arbol%centro_masa(1), Arbol%centro_masa(2))

		Arbol%masa = Arbol%masa + p%masa
		
		!OJO, deberia actualizar densidad?
		Arbol%centro_masa(0) = tmp_centromasas(0)
		Arbol%centro_masa(1) = tmp_centromasas(1)
		Arbol%centro_masa(2) = tmp_centromasas(2)

		!Tengo que redistribuir tanto a la nueva particula p, como a la que existia -> Arbol%id_particula
		!Para eso busco en cual de los hijos recien creados está contenida
		
		!write(*,*) "Rescatando informacion de particula existente q. Id: ", Arbol%id_particula, ", densidad: ", Arbol%densidad

		q = construirParticula(Arbol%id_particula, masas(Arbol%id_particula), coordenadas_x(Arbol%id_particula), coordenadas_y(Arbol%id_particula), coordenadas_z(Arbol%id_particula), 0.0D+0, 0.0D+0, 0.0D+0, Arbol%densidad)
				
		do j = 0, 7, 1			
			!obtengo el cubo que corresponde a cada hijo
						
			Hijo => DarHijo(ArbolPtr, j)

			tmp_cube = CrearCubo(Hijo%radio, Hijo%centroide)
												
			!Si el cubo hijo contiene a la particula que quiero insertar
				
			if(.not. encontrado_nuevo) then
				if(CuboContienePunto(tmp_cube, p%posicion)) then
					call InsertarParticula(p, Hijo, NNodos, masas, coordenadas_x, coordenadas_y, coordenadas_z, N, NodosParticulas)
					encontrado_nuevo = .true.
				endif
			endif
			
			if(.not. encontrado_existente) then
				if(CuboContienePunto(tmp_cube, q%posicion)) then				
					!Con esto, indico que es el nodo actual -i-, se convierte en un internal node
					Arbol%id_particula = -1
					Arbol%hoja = .false.
					Arbol%densidad = 0.0D+0
					call InsertarParticula(q, Hijo, NNodos, masas, coordenadas_x, coordenadas_y, coordenadas_z, N, NodosParticulas)
					encontrado_existente = .true.
				endif
			endif
			
			if(encontrado_existente .and. encontrado_nuevo) then
				return
			endif
		enddo		
		write(*,*) "OJO, en nodo externo no se encontro la nueva o la existente: ", encontrado_nuevo, encontrado_existente, "P y Q id's: ", p%id, q%id, "ID nodo:", Arbol%id, "Posicion particula nueva: ", p%posicion, "Posicion de particula existente: ", q%posicion
		!Debugging de problema de no deteccion en 10/13/2011
		if(.false.) then

			write(*,*) "--------------------------"
			write(*,*) "Cubo de nodo actual"
			call imprimirCubo(cube)
			write(*,*) "---------------------------------"
			write(*,*) "Info de p: (id, posicion) ", p%id, p%posicion
			write(*,*) "Info de q: (id, posicion) ", q%id, q%posicion
			write(*,*) "Contiene a p? ", CuboContienePunto(cube, p%posicion), "Radio y centroide: ", Arbol%radio, Arbol%centroide
			write(*,*) "Contiene a q? ", CuboContienePunto(cube, q%posicion)
			write(*,*) "Info de hijos"
			do j = 0, 7, 1
				write(*,*) "Info de hijo: ", Hijo%id, Hijo%id_particula, Hijo%n_particulas
				tmp_cube = CrearCubo(Hijo%radio, Hijo%centroide)
				write(*,*) "Hijo Contiene a p? ", CuboContienePunto(tmp_cube, p%posicion), "Radio y centroide: ", Hijo%radio, Hijo%centroide
				write(*,*) "Hijo Contiene a q? ", CuboContienePunto(tmp_cube, q%posicion)
			enddo			

		endif
		write(*,*) "MUERE"				
		call MPI_abort (1)

		j = 0
		!mato al programa para analizar
		j = 10 / j			
		
	endif

END SUBROUTINE

!Calcula las coordenadas tridimensionales del cubo que contiene a una esfera de radio R

function CrearCubo(radio, centroide) result(cube)
	
	use Tipos
	
	implicit none
	
	real*8 centroide(0:2)
	real*8 radio
	type(Cubo) cube
	
	cube%radio = radio
	
	cube%centroide(0) = centroide(0)
	cube%centroide(1) = centroide(1)
	cube%centroide(2) = centroide(2)
	
	cube%vertice1(0) = (-1.0 * radio) + centroide(0)
	cube%vertice1(1) = (-1.0 * radio) + centroide(1)
	cube%vertice1(2) = (radio) + centroide(2)

	cube%vertice2(0) = (-1.0 * radio) + centroide(0)
	cube%vertice2(1) = (radio) + centroide(1)
	cube%vertice2(2) = (radio) + centroide(2)

	cube%vertice3(0) = (radio) + centroide(0)
	cube%vertice3(1) = (radio) + centroide(1)
	cube%vertice3(2) = (radio) + centroide(2)

	cube%vertice4(0) = (radio) + centroide(0)
	cube%vertice4(1) = (-1.0 * radio) + centroide(1)
	cube%vertice4(2) = (radio) + centroide(2)
	
	cube%vertice5(0) = (-1.0 * radio) + centroide(0)
	cube%vertice5(1) = (-1.0 * radio) + centroide(1)
	cube%vertice5(2) = (-1.0 * radio) + centroide(2)
	
	cube%vertice6(0) = (-1.0 * radio) + centroide(0)
	cube%vertice6(1) = (radio) + centroide(1)
	cube%vertice6(2) = (-1.0 * radio) + centroide(2)
	
	cube%vertice7(0) = (radio) + centroide(0)
	cube%vertice7(1) = (radio) + centroide(1)
	cube%vertice7(2) = (-1.0 * radio) + centroide(2)

	cube%vertice8(0) = (radio) + centroide(0)
	cube%vertice8(1) = (-1.0 * radio) + centroide(1)
	cube%vertice8(2) = (-1.0 * radio) + centroide(2)
	
	return
	
end function 

!Divide el cubo de radio en centroide en 8 hijos. Requerido para el algoritmo
function DividirCubo(padre) result(hijos)
	
	use Tipos
	implicit none
	
	real*8 radio
	type(Cubo) hijos(0:7)
	type(Cubo) padre
	
	real*8 tmp_centroide(0:2)

	radio = padre%radio

	tmp_centroide(0) = padre%centroide(0) + radio/2.0
	tmp_centroide(1) = padre%centroide(1) + radio/2.0
	tmp_centroide(2) = padre%centroide(2) + radio/2.0
	
	hijos(0) = CrearCubo(radio / 2.0, tmp_centroide)
			
	tmp_centroide(0) = padre%centroide(0) + (-1.0 * radio/2.0)
	tmp_centroide(1) = padre%centroide(1) + radio/2.0
	tmp_centroide(2) = padre%centroide(2) + radio/2.0
	
	hijos(1) = CrearCubo(radio / 2.0, tmp_centroide)

	tmp_centroide(0) = padre%centroide(0) + radio/2.0
	tmp_centroide(1) = padre%centroide(1) + (-1.0 * radio/2.0)
	tmp_centroide(2) = padre%centroide(2) + radio/2.0
	
	hijos(2) = CrearCubo(radio / 2.0, tmp_centroide)

	tmp_centroide(0) = padre%centroide(0) + radio/2.0
	tmp_centroide(1) = padre%centroide(1) + radio/2.0
	tmp_centroide(2) = padre%centroide(2) + (-1.0 * radio/2.0)
	
	hijos(3) = CrearCubo(radio / 2.0, tmp_centroide)



	tmp_centroide(0) = padre%centroide(0) + (-1.0 * radio/2.0)
	tmp_centroide(1) = padre%centroide(1) + (-1.0 * radio/2.0)
	tmp_centroide(2) = padre%centroide(2) + radio/2.0
	
	hijos(4) = CrearCubo(radio / 2.0, tmp_centroide)
			
	tmp_centroide(0) = padre%centroide(0) + radio/2.0
	tmp_centroide(1) = padre%centroide(1) + (-1.0 * radio/2.0)
	tmp_centroide(2) = padre%centroide(2) + (-1.0 * radio/2.0)
	
	hijos(5) = CrearCubo(radio / 2.0, tmp_centroide)

	tmp_centroide(0) = padre%centroide(0) + (-1.0 * radio/2.0)
	tmp_centroide(1) = padre%centroide(1) + radio/2.0
	tmp_centroide(2) = padre%centroide(2) + (-1.0 * radio/2.0)
	
	hijos(6) = CrearCubo(radio / 2.0, tmp_centroide)

	tmp_centroide(0) = padre%centroide(0) + (-1.0 * radio/2.0)
	tmp_centroide(1) = padre%centroide(1) + (-1.0 * radio/2.0)
	tmp_centroide(2) = padre%centroide(2) + (-1.0 * radio/2.0)
	
	hijos(7) = CrearCubo(radio / 2.0, tmp_centroide)	
	
	return

end function 



!Busco los puntos minimos y maximos de cada eje entre los 8 vertices, si cada coordenada del punto esta en medio 
!de ambos para los 3 ejes, está dentro
!PASA A MODULO DE GEOMETRIA O FUNCIONES
logical function CuboContienePunto(cube, punto)
	use Tipos
	implicit none
	type(Cubo) cube
	real*8 punto(0:2), vertices(0:7)
	real*8 coord_min, coord_max
	logical dentro_x, dentro_y, dentro_z
	
	dentro_x = .false. 
	dentro_y = .false. 
	dentro_z = .false.
	
	!---Eje X
	vertices(0) = cube%vertice1(0)
	vertices(1) = cube%vertice2(0)
	vertices(2) = cube%vertice3(0)
	vertices(3) = cube%vertice4(0)
	vertices(4) = cube%vertice5(0)
	vertices(5) = cube%vertice6(0)
	vertices(6) = cube%vertice7(0)
	vertices(7) = cube%vertice8(0)
	
	coord_min = minval(vertices)
	coord_max = maxval(vertices)
	
	if(punto(0) <= coord_max .and.  punto(0) >= coord_min) then
		dentro_x = .true.
	endif
	
	!---Eje Y
	vertices(0) = cube%vertice1(1)
	vertices(1) = cube%vertice2(1)
	vertices(2) = cube%vertice3(1)
	vertices(3) = cube%vertice4(1)
	vertices(4) = cube%vertice5(1)
	vertices(5) = cube%vertice6(1)
	vertices(6) = cube%vertice7(1)
	vertices(7) = cube%vertice8(1)
	
	coord_min = minval(vertices)
	coord_max = maxval(vertices)
	
	if (punto(1) <= coord_max .and. punto(1) >= coord_min) then
		dentro_y = .true.
	endif

	!---Eje Y
	vertices(0) = cube%vertice1(2)
	vertices(1) = cube%vertice2(2)
	vertices(2) = cube%vertice3(2)
	vertices(3) = cube%vertice4(2)
	vertices(4) = cube%vertice5(2)
	vertices(5) = cube%vertice6(2)
	vertices(6) = cube%vertice7(2)
	vertices(7) = cube%vertice8(2)
	
	coord_min = minval(vertices)
	coord_max = maxval(vertices)
	
	if(punto(2) <= coord_max .and. punto(2) >= coord_min) then
		dentro_z = .true.
	endif
	
	!-----
	
	CuboContienePunto = dentro_x .and. dentro_y .and. dentro_z
	
	return
	
end function 



RECURSIVE SUBROUTINE ExpandirNodo(Nodo, NodoAExplorar, N, detectados, n_vecinos, lista_vecinos, distancias_vecinos)
	use Tipos
	use Constantes
	use Vectores
	use Fisica
	use Auxiliar
	
	implicit none

	integer*4 N, n_vecinos, detectados, i, j
	integer*4 lista_vecinos(0:n_vecinos-1)
	real*8 distancias_vecinos(0:n_vecinos-1)
	type(OctreeNode), TARGET :: Nodo, NodoAExplorar
	type(OctreeNode), POINTER :: dummy, NodoAExplorarPtr
	
	NodoAExplorarPtr => NodoAExplorar
	
	do j = 0, 7, 1
		dummy => DarHijo(NodoAExplorarPtr, j)
		if(associated(dummy)) then
			if(dummy%hoja .and. dummy%id_particula >= 0 .and. .not. dummy%visitado .and. .not. dummy%id == Nodo%id) then
				if(detectados < n_vecinos) then
					if(.not. dummy%visitado) then
						lista_vecinos(detectados) = dummy%id_particula
						distancias_vecinos(detectados) = distanciaPuntos(Nodo%centro_masa, dummy%centro_masa)
						detectados = detectados + 1
						dummy%visitado = .true.
					endif
				else
					!write(*,*) "RevisarReemplazar (2) en hijo: ", j
					call RevisarReemplazar(Nodo, dummy ,n_vecinos, lista_vecinos, distancias_vecinos)
				endif
			else
				if(.not. dummy%hoja .and. .not. dummy%visitado .and. dummy%n_particulas > 0) then
					call ExpandirNodo(Nodo, dummy, N, detectados, n_vecinos, lista_vecinos, distancias_vecinos)
				endif
			endif
		endif
	enddo	
	
END SUBROUTINE

SUBROUTINE RevisarReemplazar(Nodo, NodoAnalizado ,n_vecinos, lista_vecinos, distancias_vecinos)

	use Tipos
	use Auxiliar
	use Vectores
	
	implicit none

	integer*4 n_vecinos, i, j, max_idx
	integer*4 lista_vecinos(0:n_vecinos-1)
	real*8 distancias_vecinos(0:n_vecinos - 1)
	real*8 distancia_max, distancia_nodo
	type(OctreeNode) Nodo, NodoAnalizado
	
	if(Nodo%id == NodoAnalizado%id .or. NodoAnalizado%visitado) then
		return
	endif
	
	distancia_max = valorMaximoIdx(distancias_vecinos, n_vecinos, max_idx)
	distancia_nodo = distanciaPuntos(Nodo%centro_masa, NodoAnalizado%centro_masa)
	
	if(distancia_nodo < distancia_max) then
		!write(*,*) "Distancia maxima: ", distancia_max, " - particula: ", lista_vecinos(max_idx), " es mayor a distancia de nodo: ", distancia_nodo, " ID: ", NodoAnalizado%id_particula, " idx distancia maxima: ", max_idx, " Ids, centros: ", Nodo%id_particula, NodoAnalizado%id_particula, Nodo%centro_masa, NodoAnalizado%centro_masa
		lista_vecinos(max_idx) = NodoAnalizado%id_particula
		distancias_vecinos(max_idx) = distancia_nodo
		NodoAnalizado%visitado = .true.
	endif
	
END SUBROUTINE

RECURSIVE SUBROUTINE Vecinos(Nodo, NodosParticulas, Arbol, N, detectados, n_vecinos, lista_vecinos, distancias_vecinos, nivel_actual, nivel_maximo, myid)

	use Tipos
	use Constantes
	use Vectores
	use Fisica
	use Auxiliar
	
	implicit none
	
	integer*4 N, n_vecinos, detectados, i, j, nivel_actual, nivel_maximo, myid
	
	logical caso_especial

	type(OctreeNode) NodosParticulas(0:N-1)
	
	integer*4 lista_vecinos(0:n_vecinos-1)
	
	integer*4 indices_vecinos(0:n_vecinos - 1)

	real*8 distancias_vecinos(0:n_vecinos - 1)
		
	type(OctreeNode), TARGET :: Nodo, Arbol
	
	type(OctreeNode), POINTER :: dummy
	
	!Reviso adyacentes buscando partículas vecinas
	
	caso_especial = .true.

	if(Arbol%visitado .or. nivel_actual >= nivel_maximo) then	
		return
	endif
	
	Arbol%visitado = .true.
	
!!	write(*,*) "Revisando vecinos de: ", Nodo%id, Nodo%id_particula, " rama actual: ", Arbol%id,  detectados, " padre: ", Nodo%padre%id, " niveles: ", nivel_actual, nivel_maximo
	
	!Orden de revision/detección (hermanos se refiere a los de Arbol, no a los de Nodo, inicialmente, Arbol y nodo son el mismo):
		!Hermanos hojas
		!Vecinos hojas
		!Expansion de hermanos y vecinos ramas
		!Si no se han completado, se sube un nivel (llamada recursiva con Padre como Arbol).
	
	!Reviso hermanos de Arbol
	
	!No reviso hermanos si se está en la raiz
	if(Arbol%nivel > 0) then	
		do j = 0, 7, 1
			dummy => DarHijo(Arbol%padre,j)
			if(associated(dummy)) then
				if(dummy%hoja .and. dummy%id_particula >= 0 .and. .not. dummy%id == Nodo%id) then
					if(detectados < n_vecinos) then
						if(.not. dummy%visitado) then
							lista_vecinos(detectados) = dummy%id_particula
							distancias_vecinos(detectados) = distanciaPuntos(Nodo%centro_masa, dummy%centro_masa)
							detectados = detectados + 1
							dummy%visitado = .true.
						endif
					else
						call RevisarReemplazar(Nodo, dummy ,n_vecinos, lista_vecinos, distancias_vecinos)
					endif
				endif
			endif
		enddo
	endif
	
	!Reviso vecinos de Arbol
	
	if(associated(Arbol%adyacente_superior)) then
		if(Arbol%adyacente_superior%hoja .and. Arbol%adyacente_superior%id_particula >= 0 .and. .not. Arbol%adyacente_superior%id_particula == Nodo%id_particula .and. .not. Arbol%adyacente_superior%visitado) then
			if(detectados < n_vecinos) then
				if(.not. Arbol%adyacente_superior%visitado) then
					lista_vecinos(detectados) = Arbol%adyacente_superior%id_particula
					distancias_vecinos(detectados) = distanciaPuntos(Nodo%centro_masa, Arbol%adyacente_superior%centro_masa)
	!!				write(*,*) "Agregando adyacente_superior como vecino (id, distancia, detectados): ", Arbol%adyacente_superior%id_particula, distancias_vecinos(detectados), detectados
					detectados = detectados + 1
					Arbol%adyacente_superior%visitado = .true.
				endif
			else
				!write(*,*) "RevisarReemplazar adyacente_superior"
				call RevisarReemplazar(Nodo, Arbol%adyacente_superior ,n_vecinos, lista_vecinos, distancias_vecinos)
			endif
		endif
	endif

	if(associated(Arbol%adyacente_inferior)) then
		if(Arbol%adyacente_inferior%hoja .and. Arbol%adyacente_inferior%id_particula >= 0 .and. .not. Arbol%adyacente_inferior%id_particula == Nodo%id_particula .and. .not. Arbol%adyacente_inferior%visitado) then
			if(detectados < n_vecinos) then
				if(.not. Arbol%adyacente_inferior%visitado) then
					lista_vecinos(detectados) = Arbol%adyacente_inferior%id_particula
					distancias_vecinos(detectados) = distanciaPuntos(Nodo%centro_masa, Arbol%adyacente_inferior%centro_masa)
	!!				write(*,*) "Agregando adyacente_inferior como vecino (id, distancia, detectados): ", Arbol%adyacente_inferior%id_particula, distancias_vecinos(detectados), detectados
					detectados = detectados + 1
					Arbol%adyacente_inferior%visitado = .true.
				endif
			else
				!write(*,*) "RevisarReemplazar adyacente_inferior"
				call RevisarReemplazar(Nodo, Arbol%adyacente_inferior ,n_vecinos, lista_vecinos, distancias_vecinos)
			endif
		endif
	endif

	if(associated(Arbol%adyacente_izquierdo)) then
		if(Arbol%adyacente_izquierdo%hoja .and. Arbol%adyacente_izquierdo%id_particula >= 0 .and. .not. Arbol%adyacente_izquierdo%id_particula == Nodo%id_particula .and. .not. Arbol%adyacente_izquierdo%visitado) then
			if(detectados < n_vecinos) then
				if(.not. Arbol%adyacente_izquierdo%visitado) then
					lista_vecinos(detectados) = Arbol%adyacente_izquierdo%id_particula
					distancias_vecinos(detectados) = distanciaPuntos(Nodo%centro_masa, Arbol%adyacente_izquierdo%centro_masa)
	!!				write(*,*) "Agregando adyacente_izquierdo como vecino (id, distancia, detectados): ", Arbol%adyacente_izquierdo%id_particula, distancias_vecinos(detectados), detectados
					detectados = detectados + 1
					Arbol%adyacente_izquierdo%visitado = .true.
				endif
			else
				!write(*,*) "RevisarReemplazar adyacente_izquierdo"
				call RevisarReemplazar(Nodo, Arbol%adyacente_izquierdo ,n_vecinos, lista_vecinos, distancias_vecinos)
			endif
		endif		
	endif

	if(associated(Arbol%adyacente_derecho)) then
		if(Arbol%adyacente_derecho%hoja .and. Arbol%adyacente_derecho%id_particula >= 0 .and. .not. Arbol%adyacente_derecho%id_particula == Nodo%id_particula .and. .not. Arbol%adyacente_derecho%visitado) then
			if(detectados < n_vecinos) then
				if(.not. Arbol%adyacente_derecho%visitado) then
					lista_vecinos(detectados) = Arbol%adyacente_derecho%id_particula
					distancias_vecinos(detectados) = distanciaPuntos(Nodo%centro_masa, Arbol%adyacente_derecho%centro_masa)
	!!				write(*,*) "Agregando adyacente_derecho como vecino (id, distancia, detectados): ", Arbol%adyacente_derecho%id_particula, distancias_vecinos(detectados), detectados
					detectados = detectados + 1
					Arbol%adyacente_derecho%visitado = .true.
				endif
			else
				!write(*,*) "RevisarReemplazar adyacente_derecho"
				call RevisarReemplazar(Nodo, Arbol%adyacente_derecho ,n_vecinos, lista_vecinos, distancias_vecinos)
			endif
		endif		
	endif

	if(associated(Arbol%adyacente_frente)) then
		if(Arbol%adyacente_frente%hoja .and. Arbol%adyacente_frente%id_particula >= 0 .and. .not. Arbol%adyacente_frente%id_particula == Nodo%id_particula .and. .not. Arbol%adyacente_frente%visitado) then
			if(detectados < n_vecinos) then
				if(.not. Arbol%adyacente_frente%visitado) then
					lista_vecinos(detectados) = Arbol%adyacente_frente%id_particula
					distancias_vecinos(detectados) = distanciaPuntos(Nodo%centro_masa, Arbol%adyacente_frente%centro_masa)
	!!				write(*,*) "Agregando adyacente_frente como vecino (id, distancia, detectados): ", Arbol%adyacente_frente%id_particula, distancias_vecinos(detectados), detectados
					detectados = detectados + 1
					Arbol%adyacente_frente%visitado = .true.
				endif
			else
				!write(*,*) "RevisarReemplazar adyacente_frente"
				call RevisarReemplazar(Nodo, Arbol%adyacente_frente ,n_vecinos, lista_vecinos, distancias_vecinos)
			endif
		endif		
	endif

	if(associated(Arbol%adyacente_trasero)) then
		if(Arbol%adyacente_trasero%hoja .and. Arbol%adyacente_trasero%id_particula >= 0 .and. .not. Arbol%adyacente_trasero%id_particula == Nodo%id_particula .and. .not. Arbol%adyacente_trasero%visitado) then
			if(detectados < n_vecinos) then
				if(.not. Arbol%adyacente_trasero%visitado) then
					lista_vecinos(detectados) = Arbol%adyacente_trasero%id_particula
					distancias_vecinos(detectados) = distanciaPuntos(Nodo%centro_masa, Arbol%adyacente_trasero%centro_masa)
	!!				write(*,*) "Agregando adyacente_trasero como vecino (id, distancia, detectados): ", Arbol%adyacente_trasero%id_particula, distancias_vecinos(detectados), detectados
					detectados = detectados + 1
					Arbol%adyacente_trasero%visitado = .true.
				endif
			else
				!write(*,*) "RevisarReemplazar adyacente_trasero"
				call RevisarReemplazar(Nodo, Arbol%adyacente_trasero ,n_vecinos, lista_vecinos, distancias_vecinos)
			endif		
		endif		
	endif		

	!------------------
	!Expando hermanos rama
	!------------------
	
!!	write(*,*) "Antes de expandir hermano rama. Detectados: ", detectados

	if(Arbol%nivel > 0) then	
		do j = 0, 7, 1
			dummy => DarHijo(Arbol%padre, j)
			if(associated(dummy)) then
				if(.not. dummy%hoja .and. dummy%n_particulas > 0) then
					call ExpandirNodo(Nodo, dummy, N, detectados, n_vecinos, lista_vecinos, distancias_vecinos)
				endif
			endif
		enddo
	endif
	
	!------------------
	!Expando adyacentes
	!------------------
			
	if(associated(Arbol%adyacente_superior)) then		
				
		if(.not. Arbol%adyacente_superior%hoja .and. Arbol%adyacente_superior%n_particulas > 0 .and. .not. Arbol%adyacente_superior%visitado) then
			call ExpandirNodo(Nodo, Arbol%adyacente_superior, N, detectados, n_vecinos, lista_vecinos, distancias_vecinos)
		endif
	endif

	if(associated(Arbol%adyacente_inferior)) then
				
		if(.not. Arbol%adyacente_inferior%hoja .and. Arbol%adyacente_inferior%n_particulas > 0 .and. .not. Arbol%adyacente_inferior%visitado) then
			call ExpandirNodo(Nodo, Arbol%adyacente_inferior, N, detectados, n_vecinos, lista_vecinos, distancias_vecinos)
		endif
	endif

	if(associated(Arbol%adyacente_izquierdo)) then
				
		if(.not. Arbol%adyacente_izquierdo%hoja .and. Arbol%adyacente_izquierdo%n_particulas > 0 .and. .not. Arbol%adyacente_izquierdo%visitado) then
			call ExpandirNodo(Nodo, Arbol%adyacente_izquierdo, N, detectados, n_vecinos, lista_vecinos, distancias_vecinos)
		endif
	endif

	if(associated(Arbol%adyacente_derecho)) then
				
		if(.not. Arbol%adyacente_derecho%hoja .and. Arbol%adyacente_derecho%n_particulas > 0 .and. .not. Arbol%adyacente_derecho%visitado) then
			call ExpandirNodo(Nodo, Arbol%adyacente_derecho, N, detectados, n_vecinos, lista_vecinos, distancias_vecinos)
		endif
	endif

	if(associated(Arbol%adyacente_frente)) then
				
		if(.not. Arbol%adyacente_frente%hoja .and. Arbol%adyacente_frente%n_particulas > 0 .and. .not. Arbol%adyacente_frente%visitado) then
			call ExpandirNodo(Nodo, Arbol%adyacente_frente, N, detectados, n_vecinos, lista_vecinos, distancias_vecinos)
		endif
	endif

	if(associated(Arbol%adyacente_trasero)) then
				
		if(.not. Arbol%adyacente_trasero%hoja .and. Arbol%adyacente_trasero%n_particulas > 0 .and. .not. Arbol%adyacente_trasero%visitado) then
			call ExpandirNodo(Nodo, Arbol%adyacente_trasero, N, detectados, n_vecinos, lista_vecinos, distancias_vecinos)
		endif
	endif
	
	!------------------
	!Llamada recursiva
	!------------------

	if(detectados < n_vecinos .or. nivel_actual < nivel_maximo) then
		
		if(associated(Arbol%adyacente_superior)) then

			if(.not. Arbol%adyacente_superior%hoja .and. Arbol%adyacente_superior%n_particulas > 0) then
				caso_especial = .false.
				call Vecinos(Nodo, NodosParticulas, Arbol%adyacente_superior, N, detectados, n_vecinos, lista_vecinos, distancias_vecinos, nivel_actual + 1, nivel_maximo, myid)
			endif
		endif

		if(associated(Arbol%adyacente_inferior)) then

			if(.not. Arbol%adyacente_inferior%hoja .and. Arbol%adyacente_inferior%n_particulas > 0) then
				caso_especial = .false.
				call Vecinos(Nodo, NodosParticulas, Arbol%adyacente_inferior, N, detectados, n_vecinos, lista_vecinos, distancias_vecinos, nivel_actual + 1, nivel_maximo, myid)
			endif
		endif

		if(associated(Arbol%adyacente_izquierdo)) then

			if(.not. Arbol%adyacente_izquierdo%hoja .and. Arbol%adyacente_izquierdo%n_particulas > 0) then
				caso_especial = .false.
				call Vecinos(Nodo, NodosParticulas, Arbol%adyacente_izquierdo, N, detectados, n_vecinos, lista_vecinos, distancias_vecinos, nivel_actual + 1, nivel_maximo, myid)
			endif
		endif

		if(associated(Arbol%adyacente_derecho)) then

			if(.not. Arbol%adyacente_derecho%hoja .and. Arbol%adyacente_derecho%n_particulas > 0) then
				caso_especial = .false.
				call Vecinos(Nodo, NodosParticulas, Arbol%adyacente_derecho, N, detectados, n_vecinos, lista_vecinos, distancias_vecinos, nivel_actual + 1, nivel_maximo, myid)
			endif
		endif

		if(associated(Arbol%adyacente_frente)) then

			if(.not. Arbol%adyacente_frente%hoja .and. Arbol%adyacente_frente%n_particulas > 0) then
				caso_especial = .false.
				call Vecinos(Nodo, NodosParticulas, Arbol%adyacente_frente, N, detectados, n_vecinos, lista_vecinos, distancias_vecinos, nivel_actual + 1, nivel_maximo, myid)
			endif
		endif

		if(associated(Arbol%adyacente_trasero)) then

			if(.not. Arbol%adyacente_trasero%hoja .and. Arbol%adyacente_trasero%n_particulas > 0) then
				caso_especial = .false.
				call Vecinos(Nodo, NodosParticulas, Arbol%adyacente_trasero, N, detectados, n_vecinos, lista_vecinos, distancias_vecinos, nivel_actual + 1, nivel_maximo, myid)
			endif
		endif		
		
		!Llamada al padre, no incrementa nivel actual en caso especial de que todos los vecinos y hermanos sean hojas
		
		if( (caso_especial .or. detectados < n_vecinos) .and. Arbol%id > 0) then
		!if(caso_especial .or. detectados < n_vecinos) then
			call Vecinos(Nodo, NodosParticulas, Arbol%padre, N, detectados, n_vecinos, lista_vecinos, distancias_vecinos, nivel_actual, nivel_maximo, myid)
		endif
	endif

END SUBROUTINE

SUBROUTINE AsignarAdyacentes(posicion_hijo, Nodo, Padre)

	use Tipos
	
	implicit none

	integer*4 posicion_hijo
	
	type(OctreeNode), TARGET :: Nodo, Padre

	type(OctreeNode), POINTER :: Hijo, PadrePtr

	PadrePtr => Padre

	SELECT CASE (posicion_hijo)
	   
	   CASE (0)	   	
	   	Nodo%adyacente_izquierdo => DarHijo(PadrePtr, 1)
	   	Nodo%adyacente_frente => DarHijo(PadrePtr, 2)
	   	Nodo%adyacente_inferior => DarHijo(PadrePtr, 3)
	   	
	   	!!Vecinos de padre
	   	
	   	if(associated(Padre%adyacente_derecho)) then
	   		Nodo%adyacente_derecho => Padre%adyacente_derecho
	   	endif	   	

	   	if(associated(Padre%adyacente_trasero)) then
	   		Nodo%adyacente_trasero => Padre%adyacente_trasero
	   	endif
	   	
	   	if(associated(Padre%adyacente_superior)) then
	   		Nodo%adyacente_superior => Padre%adyacente_superior
	   	endif
	   	
	   CASE (1)		
	   	Nodo%adyacente_derecho => DarHijo(PadrePtr, 0)
	   	Nodo%adyacente_frente => DarHijo(PadrePtr, 4)
	   	Nodo%adyacente_inferior => DarHijo(PadrePtr, 6)
	   	
	   	!!Vecinos de padre
	   	if(associated(Padre%adyacente_izquierdo)) then
	   		Nodo%adyacente_izquierdo => Padre%adyacente_izquierdo
	   	endif
	   	
	   	if(associated(Padre%adyacente_trasero)) then
	   		Nodo%adyacente_trasero => Padre%adyacente_trasero
	   	endif
	   	
	   	if(associated(Padre%adyacente_superior)) then
	   		Nodo%adyacente_superior => Padre%adyacente_superior
	   	endif
	   	
	   	
	   CASE (2)
	   	Nodo%adyacente_izquierdo => DarHijo(PadrePtr, 4)
	   	Nodo%adyacente_trasero => DarHijo(PadrePtr, 0)
	   	Nodo%adyacente_inferior => DarHijo(PadrePtr, 5)
	   	
	   	!!Vecinos de padre
	   	if(associated(Padre%adyacente_derecho)) then
	   		Nodo%adyacente_derecho => Padre%adyacente_derecho
	   	endif
	   	
	   	if(associated(Padre%adyacente_frente)) then
	   		Nodo%adyacente_frente => Padre%adyacente_frente
	   	endif
	   	
	   	if(associated(Padre%adyacente_superior)) then
	   		Nodo%adyacente_superior => Padre%adyacente_superior
	   	endif
	   		   	
	   CASE (3)
	   	Nodo%adyacente_izquierdo => DarHijo(PadrePtr, 6)
	   	Nodo%adyacente_frente => DarHijo(PadrePtr, 5)
	   	Nodo%adyacente_superior => DarHijo(PadrePtr, 0)
	   	
	   	!!Vecinos de padre
	   	if(associated(Padre%adyacente_derecho)) then
	   		Nodo%adyacente_derecho => Padre%adyacente_derecho
	   	endif
	   	
	   	if(associated(Padre%adyacente_trasero)) then
	   		Nodo%adyacente_trasero => Padre%adyacente_trasero
	   	endif
	   	
	   	if(associated(Padre%adyacente_inferior)) then
	   		Nodo%adyacente_inferior => Padre%adyacente_inferior
	   	endif	   	
	   	
	   CASE (4)
	   	Nodo%adyacente_derecho => DarHijo(PadrePtr, 2)
	   	Nodo%adyacente_trasero => DarHijo(PadrePtr, 1)
	   	Nodo%adyacente_inferior => DarHijo(PadrePtr, 7)	   
	   	
	   	!!Vecinos de padre
	   	if(associated(Padre%adyacente_izquierdo)) then
	   		Nodo%adyacente_izquierdo => Padre%adyacente_izquierdo
	   	endif
	   	
	   	if(associated(Padre%adyacente_frente)) then
	   		Nodo%adyacente_frente => Padre%adyacente_frente
	   	endif
	   	
	   	if(associated(Padre%adyacente_superior)) then
	   		Nodo%adyacente_superior => Padre%adyacente_superior
	   	endif	   	
	   	
	   	
	   CASE (5)
	   	Nodo%adyacente_izquierdo => DarHijo(PadrePtr, 7)
	   	Nodo%adyacente_trasero => DarHijo(PadrePtr, 3)
	   	Nodo%adyacente_superior => DarHijo(PadrePtr, 2)
	   	
	   	!!Vecinos de padre
	   	if(associated(Padre%adyacente_derecho)) then
	   		Nodo%adyacente_derecho => Padre%adyacente_derecho
	   	endif
	   	
	   	if(associated(Padre%adyacente_frente)) then
	   		Nodo%adyacente_frente => Padre%adyacente_frente
	   	endif
	   	
	   	if(associated(Padre%adyacente_inferior)) then
	   		Nodo%adyacente_inferior => Padre%adyacente_inferior
	   	endif	   	
	   	
	   	
	   CASE (6)
	   	Nodo%adyacente_derecho => DarHijo(PadrePtr, 3)
	   	Nodo%adyacente_frente => DarHijo(PadrePtr, 7)
	   	Nodo%adyacente_superior => DarHijo(PadrePtr, 1)
	   	
	   	!!Vecinos de padre
	   	if(associated(Padre%adyacente_izquierdo)) then
	   		Nodo%adyacente_izquierdo => Padre%adyacente_izquierdo
	   	endif
	   	
	   	if(associated(Padre%adyacente_trasero)) then
	   		Nodo%adyacente_trasero => Padre%adyacente_trasero
	   	endif
	   	
	   	if(associated(Padre%adyacente_inferior)) then
	   		Nodo%adyacente_inferior => Padre%adyacente_inferior
	   	endif	   	
	   	
	   	
	   CASE (7)
	   	Nodo%adyacente_derecho => DarHijo(PadrePtr, 5)
	   	Nodo%adyacente_trasero => DarHijo(PadrePtr, 6)
	   	Nodo%adyacente_superior => DarHijo(PadrePtr, 4)	   	   
	      
	   	if(associated(Padre%adyacente_izquierdo)) then
	   		Nodo%adyacente_izquierdo => Padre%adyacente_izquierdo
	   	endif
	   	
	   	if(associated(Padre%adyacente_frente)) then
	   		Nodo%adyacente_frente => Padre%adyacente_frente
	   	endif
	   	
	   	if(associated(Padre%adyacente_inferior)) then
	   		Nodo%adyacente_inferior => Padre%adyacente_inferior
	   	endif	   	
	      
	      
	END SELECT

END SUBROUTINE

RECURSIVE SUBROUTINE VecinosFuerzaBruta(NodosParticulas, IdNodo, N, n_vecinos, lista_vecinos)

	use Funciones
	use Auxiliar
	use Tipos
	use Vectores

	implicit none

	type(OctreeNode) NodosParticulas(0:N-1)
	integer*4 lista_vecinos(n_vecinos-1)
	real*8 distancias(0:N-1)
	integer*4 i, N, n_vecinos, IdNodo
	integer*4 indices_vecinos(0:N-1)
	
	do i = 0, N - 1, 1
		indices_vecinos(i) = i
		distancias(i) = distanciaPuntos(NodosParticulas(IdNodo)%centro_masa, NodosParticulas(i)%centro_masa)
	enddo

	distancias = ordenar(distancias, indices_vecinos, N)
	
	!Empiezo del elemento 1 porque el elemento cero será el mismo (distancia cero)
!	do i = 1, n_vecinos, 1
!		lista_vecinos(i - 1) = NodosParticulas(indices_vecinos(i))%id_particula
!	enddo			

	do i = 0, n_vecinos - 1, 1
		lista_vecinos(i) = NodosParticulas(indices_vecinos(i))%id_particula
	enddo			
		

END SUBROUTINE

RECURSIVE SUBROUTINE ResetearVisitados(Arbol)

	use Tipos
	
	implicit none

	integer*4 i
		
	type(OctreeNode), TARGET :: Arbol

	type(OctreeNode), POINTER :: ArbolPtr

	ArbolPtr => Arbol

	Arbol%visitado = .false.
	
	if(.not. Arbol%hoja) then
		
		do i = 0, 7, 1
			call ResetearVisitados(DarHijo(ArbolPtr, i))
		enddo
	endif
	

END SUBROUTINE

subroutine imprimirCubo(cube)
	use Tipos
	implicit none
	type(Cubo) cube
	
	write(*,*) cube%vertice1(0), ",", cube%vertice1(1), ",", cube%vertice1(2)
	write(*,*) cube%vertice2(0), ",", cube%vertice2(1), ",", cube%vertice2(2)
	write(*,*) cube%vertice3(0), ",", cube%vertice3(1), ",", cube%vertice3(2)
	write(*,*) cube%vertice4(0), ",", cube%vertice4(1), ",", cube%vertice4(2)
	write(*,*) cube%vertice5(0), ",", cube%vertice5(1), ",", cube%vertice5(2)
	write(*,*) cube%vertice6(0), ",", cube%vertice6(1), ",", cube%vertice6(2)
	write(*,*) cube%vertice7(0), ",", cube%vertice7(1), ",", cube%vertice7(2)
	write(*,*) cube%vertice8(0), ",", cube%vertice8(1), ",", cube%vertice8(2)
	
end subroutine 

function DarHijo(Arbol, N) result(hijo)
	
	use Tipos
	
	implicit none
	
	type(OctreeNode), POINTER :: Arbol, hijo
	
	integer*4 N

	SELECT CASE (N)
	   
	   CASE (0)	   	
		hijo => Arbol%hijo0

	   CASE (1)	   	
		hijo => Arbol%hijo1

	   CASE (2)
		hijo => Arbol%hijo2

	   CASE (3)	   	
		hijo => Arbol%hijo3

	   CASE (4)	   	
		hijo => Arbol%hijo4

	   CASE (5)	   	
		hijo => Arbol%hijo5

	   CASE (6)	   	
		hijo => Arbol%hijo6

	   CASE (7)	   	
		hijo => Arbol%hijo7
	
	END SELECT

	return
	
end function 


END MODULE Octree
