MODULE Vectores
	CONTAINS

!-----------------------------------------------
!Funcion auxiliar para calcular magnitud de vector 3D
!-----------------------------------------------

real*8 function magnitudVector3D(vector)
	
	implicit none
	real*8 vector(0:2)	

	magnitudVector3D = SQRT(vector(0)**2 + vector(1)**2 + vector(2)**2)

	return

end function

real*8 function magnitudVector2D(vector)
	
	implicit none
	real*8 vector(0:1)	

	magnitudVector2D = SQRT(vector(0)**2 + vector(1)**2)

	return

end function


function diferenciaVectores3D(vector1, vector2) result(diferencia)
	
	implicit none
	real*8 vector1(0:2), vector2(0:2), diferencia(0:2)

	diferencia(0) = vector1(0) - vector2(0)
	diferencia(1) = vector1(1) - vector2(1)
	diferencia(2) = vector1(2) - vector2(2)

	return

end function

real*8 function distanciaPuntos(p1, p2)
	use Constantes
	use Tipos

	implicit none

	real*8 p1(0:2), p2(0:2)

	distanciaPuntos = SQRT( ((p1(0) - p2(0))**2) + ((p1(1) - p2(1))**2) + ((p1(2) - p2(2))**2) )
end function

!Producto cruz de 2 vectores 3D
function productoCruz3D(vector1, vector2) result(producto)
	
	implicit none
	real*8 vector1(0:2), vector2(0:2), producto(0:2)

	producto(0) = (vector1(1) * vector2(2)) - (vector1(2) * vector2(1))
	producto(1) = (vector1(2) * vector2(0)) - (vector1(0) * vector2(2))
	producto(2) = (vector1(0) * vector2(1)) - (vector1(1) * vector2(0))	

	return
end function

!Producto escalar de 2 vectores 3D
real*8 function productoEscalar3D(vector1, vector2)
	
	implicit none
	real*8 vector1(0:2), vector2(0:2)

	productoEscalar3D = (vector1(0) * vector2(0)) + (vector1(1) * vector2(1)) + (vector1(2) * vector2(2))

	return
end function


END MODULE Vectores
