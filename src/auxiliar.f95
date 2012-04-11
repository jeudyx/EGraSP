MODULE Auxiliar
	CONTAINS

!Funcion para ordenar un arreglo
!arreglo: arreglo a ordenar
!N: tamaño del arreglo
!resultado: respuesta

recursive function ordenar(arreglo, indices, N) result(resultado)
	implicit none
		
	integer N, N1, N2
	
	real*8 arreglo(0:N-1), resultado(0:N-1)
	integer indices(0:N-1), resp_indices(0:N-1)

	resp_indices = 0.0D+0
	
	!write(*,*) "Hola, estoy en Auxiliar.ordenar. N = ", N, ", arreglo(0) = ", arreglo(0), ", indices(0) = ", indices(0)
	
	if(N == 1) then
		resultado(0) = arreglo(0)
		return
	else if (N == 2) then
		if(arreglo(0) <= arreglo(1)) then
			resultado(0) = arreglo(0)
			resultado(1) = arreglo(1)
			resp_indices(0) = indices(0)
			resp_indices(1) = indices(1)
		else
			resultado(1) = arreglo(0)
			resultado(0) = arreglo(1)		
			resp_indices(1) = indices(0)
			resp_indices(0) = indices(1)			
		endif
		
		indices = resp_indices
		return
	else
		
		N1 = N / 2
		N2 = (N - N1)
		resultado = appendArrayOrdenado(ordenar(arreglo(0:N1-1), indices(0:N1-1), N1), ordenar(arreglo(N1:N-1), indices(N1:N-1), N2), N1, N2, indices(0:N1-1), indices(N1:N-1), resp_indices)
		indices = resp_indices
		return
	endif
		
end function

function appendArrayOrdenado(arr1, arr2, N1, N2, idx1, idx2, indices) result(resultado)

	implicit none

	real*8 arr1(0:N1-1), arr2(0:N2-1) 
	integer idx1(0:N1-1), idx2(0:N2-1), indices(0:(N1+N2-1))
	real*8 resultado(0:(N1+N2-1))
	
	integer i, j, N1, N2, N, ptr1, ptr2
	
	logical seguir
	
	seguir = .true.

	N = N1 + N2 !max(N1, N2)
	
	ptr1 = 0
	ptr2 = 0
	i = 0
	resultado = 0

!	write(*,*), 'En append ordenado: ', idx1, '--', idx2, '--'
	
	do while(i <= N .and. seguir)
		if(ptr1 >= N1) then
!			write(*,*) '-', i, ptr1, ptr2, i + N2 - ptr2, '--', arr2(ptr2: N2 - 1), '--', N2 - ptr2, '__', resultado(0:i-1)
			seguir = .false.	
			resultado = appendArrayReal(resultado(0:i-1), arr2(ptr2: N2 - 1), i, N2 - ptr2)
			indices = appendArrayInt(indices(0:i-1), idx2(ptr2: N2 - 1), i, N2 - ptr2)
		else if(ptr2 >= N2) then
			seguir = .false.	
			resultado = appendArrayReal(resultado(0:i-1), arr1(ptr1: N1 - 1), i, N1 - ptr1)
			indices = appendArrayInt(indices(0:i-1), idx1(ptr1: N1 - 1), i, N1 - ptr1)
		else
			if(arr1(ptr1) < arr2(ptr2)) then
				resultado(i) = arr1(ptr1)
				indices(i) = idx1(ptr1)
				ptr1 = ptr1 + 1
			else
				resultado(i) = arr2(ptr2)
				indices(i) = idx2(ptr2)
				ptr2 = ptr2 + 1			
			endif
		endif
		i = i + 1
	end do
	
!	write(*,*), 'En append ordenado: ', indices
	
	return

end function


function appendArrayReal(arr1, arr2, N1, N2) result(resultado)

	implicit none

	real*8 arr1(0:N1-1), arr2(0:N2-1)
	
	real*8 resultado(0:(N1+N2-1))
	
	integer i, j, N1, N2
		
	do i = 0, N1 - 1, 1
		resultado(i) = arr1(i)
	end do

	do j = 0, N2 - 1, 1
		resultado(j+i) = arr2(j)
	end do

	return

end function

function appendArrayInt(arr1, arr2, N1, N2) result(resultado)

	implicit none

	integer arr1(0:N1-1), arr2(0:N2-1)
	
	integer resultado(0:(N1+N2-1))
	
	integer i, j, N1, N2
		
	do i = 0, N1 - 1, 1
		resultado(i) = arr1(i)
	end do

	do j = 0, N2 - 1, 1
		resultado(j+i) = arr2(j)
	end do

	return

end function

function appendArrayOctree(arr1, arr2, N1, N2) result(resultado)

	use Tipos

	implicit none
	
	type(OctreeNode) arr1(0:N1-1), arr2(0:N2-1)
	
	type(OctreeNode) resultado(0:(N1+N2-1))
	
	integer i, j, N1, N2
		
	do i = 0, N1 - 1, 1
		resultado(i) = arr1(i)
	end do

	do j = 0, N2 - 1, 1
		resultado(j+i) = arr2(j)
	end do

	return

end function

real*8 function logaritmo(base, y)

	implicit none
	
	real*8 y, base
	
	logaritmo = log(y) / log(base)

	return

end function

real*8 function valorMaximoIdx(arreglo, N, indice_mayor)
	implicit none
	
	integer N, indice_mayor, i
	real*8 tmp_max
	real*8 arreglo(0: N - 1)

	tmp_max = arreglo(0)
	indice_mayor = 0
	
	do i = 0, N - 1, 1
		if(arreglo(i) > tmp_max) then
			tmp_max = arreglo(i)
			indice_mayor = i
		endif
	enddo

	valorMaximoIdx = tmp_max

	return

end function

real*8 function valorMinimoIdx(arreglo, N, indice_menor)
	implicit none
	
	integer N, indice_menor, i
	real*8 tmp_min
	real*8 arreglo(0: N - 1)

	tmp_min = arreglo(0)
	indice_menor = 0
	
	do i = 0, N - 1, 1
		if(arreglo(i) < tmp_min) then
			tmp_min = arreglo(i)
			indice_menor = i
		endif
	enddo

	valorMinimoIdx = tmp_min

	return

end function

END MODULE Auxiliar
