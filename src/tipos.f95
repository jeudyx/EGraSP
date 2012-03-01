module Tipos

	type runparameters
		character (len=256) :: title
		character (len=64) :: institution
		character (len=64) :: source
		character (len=128) :: comment
		double precision :: dt
		double precision :: totalmass
		double precision :: initial_density
		double precision :: beta
		double precision :: temperature
		double precision :: BH_theta
		integer :: N_neighbour
		integer :: save_every
	end type runparameters
		
	type Particula	
		real*8 posicion(0:2)
		real*8 velocidad(0:2)			
		integer :: id
		real*8 :: masa
		real*8 :: densidad
	end type Particula
	
	type Capa
		integer :: id 		
		real*8 :: masa
		real*8 :: radio		
		integer :: num_particulas
	end type Capa
	
	type OctreeNode
		integer id
		integer id_particula
		integer n_particulas
		integer nivel
		logical hoja
		logical hijos_creados
		logical visitado
		type(OctreeNode), POINTER :: hijos(:)
		type(OctreeNode), POINTER :: padre
		real*8 centroide(0:2)
		real*8 centro_masa(0:2)
		real*8 radio
		real*8 masa
		real*8 densidad
		type(OctreeNode), POINTER :: adyacente_superior => NULL()
		type(OctreeNode), POINTER :: adyacente_inferior => NULL() 
		type(OctreeNode), POINTER :: adyacente_izquierdo => NULL()
		type(OctreeNode), POINTER :: adyacente_derecho => NULL()
		type(OctreeNode), POINTER :: adyacente_frente => NULL()
		type(OctreeNode), POINTER :: adyacente_trasero => NULL()
	end type OctreeNode

	type Cubo
		real*8 radio
		real*8 centroide(0:2)
		real*8 vertice1(0:2)
		real*8 vertice2(0:2)
		real*8 vertice3(0:2)
		real*8 vertice4(0:2)
		real*8 vertice5(0:2)
		real*8 vertice6(0:2)
		real*8 vertice7(0:2)
		real*8 vertice8(0:2)
	end type Cubo


end module Tipos
