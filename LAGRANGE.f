C
C       LAGRANGIAN INTERPOLATING POLYNOMIAL
C       FROM STEVEN CHAPRA AND TIM WOOL, 12/8/2020
C
    	REAL FUNCTION alagrange(x, y, order)
    	! input:
    	!   x = vector of independent variables (depth)
    	!   y = vector of dependent variables (oxygen concentration)
    	!   order = order of polynomial (linear = 1, quadratic = 2)
    	! output:
    	!   alagrange = dependent variable at x = 0 (interface oxygen)

    	INTEGER :: i, j, order
    	REAL :: sum, product
    	REAL :: x(3), y(3)
    	nn = order + 1
    	sum = 0
    	DO i = 1, nn
      	product = y(i)
      	DO j = 1, nn
        		IF (i .ne. j) THEN
          		product = product * (-x(j)) / (x(i) - x(j))
        		END IF
      	END DO
      	sum = sum + product
    	END DO
    	alagrange = sum

    	END FUNCTION alagrange