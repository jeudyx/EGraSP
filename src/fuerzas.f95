MODULE Fuerzas
	CONTAINS

!Modulo que contiene el calculo de fuerzas/aceleraciones (gravedad, presion, campo magnetico, etc)

function calcularAceleracion(p, Arbol, umbralBH, soft_len, N) result(acc_vect)

	use BarnesHut
	use Tipos
	use Octree
	
	implicit none
	
	real*8 acc_vect(0:2)
	!type(OctreeNode) Arbol
	
	type(OctreeNode), POINTER :: Arbol
	
	real*8 umbralBH, soft_len
		
	integer*4 j, i, N	
	
	type(Particula) p
	
	
	acc_vect(0) = 0.0D+0
	acc_vect(1) = 0.0D+0
	acc_vect(2) = 0.0D+0	
	
	if(p%masa == 0.0D+0 .or. p%densidad == 0.0D+0) then
		!Particulas eliminadas no tienen aceleracion
		return
	else
		!De momento tengo solo la parte gravitacional (10/7/2011)
		acc_vect = calcularAccGravBH(p, Arbol, umbralBH, soft_len)
		
		return	
	endif
		
end function

END MODULE Fuerzas


