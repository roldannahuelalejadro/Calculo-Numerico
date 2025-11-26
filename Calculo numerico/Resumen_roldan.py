import numpy as np
"""
METODOS:

metodo_euler(x0, t0, f , h , N):...

ÁLGEBRA LINEAL Y MATRICES:
np.linalg.solve(A, b)        # Resuelve sistema lineal Ax = b
np.linalg.inv(A)             # Calcula inversa de matriz A
np.linalg.eig(A)             # Calcula autovalores y autovectores
np.linalg.eigvals(A)         # Calcula solo autovalores
np.linalg.qr(A)              # Descomposición QR de matriz A
np.linalg.norm(x, ord)       # Norma vector/matriz (1, 2, inf)
np.linalg.cond(A)            # Número de condición de matriz A
np.linalg.det(A)             # Determinante de matriz A
np.linalg.matrix_rank(A)     # Rango de matriz A

MÉTODOS ITERATIVOS:
np.diag(A)                   # Extrae diagonal de matriz
np.diagflat(v)               # Crea matriz diagonal desde vector
np.tril(A), np.triu(A)       # Partes triangular inferior/superior
np.diag_indices(n)           # Índices para diagonal principal

INTERPOLACIÓN Y AJUSTE:
np.vander(x, N)              # Matriz de Vandermonde para polinomios
np.polyfit(x, y, grado)      # Ajuste polinomial por mínimos cuadrados
np.polyval(p, x)             # Evalúa polinomio con coeficientes p (usa base de monomios [a,b,c,d] =  )
np.polyder(p)                # Derivada de polinomio

INTEGRACIÓN NUMÉRICA:
np.trapz(y, x)               # Integración por regla trapezoidal
np.cumsum(y) * h             # Integración rectangular acumulada

BÚSQUEDA DE RAÍCES:
np.roots(p)                  # Raíces de polinomio con coeficientes p
np.sign(x)                   # Signo (útil para método bisección)
np.abs(x)                    # Valor absoluto (para criterio parada)

PROPIEDADES MATRICIALES:
np.all(eigvals > 0)          # Verifica matriz definida positiva
np.sum(np.abs(A), axis=1)    # Suma filas para dominancia diagonal
np.trace(A)                  # Traza de matriz (suma diagonal)

CONSTRUCCIÓN DE MATRICES:
np.eye(n)                    # Matriz identidad n×n
np.zeros((n, m))             # Matriz de ceros
np.ones((n, m))              # Matriz de unos
np.diag(v)                   # Matriz diagonal desde vector

OPERACIONES VECTORIALES:
np.dot(a, b)                 # Producto punto/punto
np.cross(a, b)               # Producto cruz (3D)
np.outer(a, b)               # Producto exterior
np.inner(a, b)               # Producto interno

GENERACIÓN DE MALLAS:
np.linspace(a, b, n)         # n puntos equiespaciados en [a,b]
np.arange(a, b, h)           # Puntos desde a hasta b con paso h
np.meshgrid(x, y)            # Malla 2D para superficies

CÁLCULO DE ERRORES:
np.abs(x_exact - x_approx)   # Error absoluto
np.abs((x_exact - x_approx)/x_exact)  # Error relativo
np.diff(y)                   # Diferencias finitas (derivadas)

MÉTODOS ESPECÍFICOS:
rho = max(abs(eigvals(B)))   # Radio espectral matriz B
h = (b-a)/n                  # Paso para integración
x[1:] - x[:-1]               # Diferencias entre puntos consecutivos

VERIFICACIÓN CONVERGENCIA:
np.allclose(A, B)                   # Compara si matrices son cercanas
np.max(np.abs(x_new - x_old))       # Norma infinito diferencia
np.linalg.matrix_power(A, n)        # Potencia de matriz (método potencia)
"""