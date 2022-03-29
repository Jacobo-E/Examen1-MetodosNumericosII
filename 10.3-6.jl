#= 
CÓDIGO DEL MÉTODO NUMÉRICO "NEWTON-RAPHSON" PARA FUNCIONES MULTIVARIABLES 

Autor: Escorcia Alcantara Gregorio Jacobo

Descripción: El código obtiene las raíces de una función en R^n.
Como usarlo: Todo el código es funcional, solo se debe modificar la 
sección de funciones, en donde puedes introducir las funciones de las 
cuales quieres obtener sus raíces. No olvidar modificar la matriz 
Jacobiana, en la cual se deben poner las derivadas correspondientes a 
las funciones.

IMPORTANTE: Si se quiere trabajar con una función vectorial de R^n con 
n diferente a 3, se deben agregar o quitar variables x_Kn y x_Km, 
además de modificar los vectores dentro de la o las funciones. No 
es necesario modificar la parte de las soluciones del sistema por el 
método de Gauss Jordan.
=#

#Constantes y Variables
const TOL = 0.000001
const max_iter = 100

# ----- FUNCIONES -----
#Funciones originales
f_1(x_1, x_2, x_3, x_4) = 4*x_1 - x_2 + x_3 - x_1*x_4
f_2(x_1, x_2, x_3, x_4) = -x_1 + 3*x_2 - 2*x_3 - x_2*x_4
f_3(x_1, x_2, x_3, x_4) = x_1 - 2*x_2 + 3*x_3 - x_3*x_4
f_3(x_1, x_2, x_3, x_4) = x_1^2 + x_2^2 + x_3^2 - 1 - 1/3

function F(x_1, x_2, x_3, x_4)

    f_1 = 4*x_1 - x_2 + x_3 - x_1*x_4
    f_2 = -x_1 + 3*x_2 - 2*x_3 - x_2*x_4
    f_3 = x_1 - 2*x_2 + 3*x_3 - x_3*x_4
    f_4 = x_1^2 + x_2^2 + x_3^2 - 1 - 1/3

    return [f_1; f_2; f_3; f_4]
end
#----------------------

#Matriz Jacobiana:
function J(x_1, x_2, x_3, x_4)
    a_11 = 4 - x_4
    a_12 = -1
    a_13 = 1
    a_14 = -x_1

    a_21 = -1
    a_22 = 3 - x_4
    a_23 = -2
    a_24 = -x_2

    a_31 = 1
    a_32 = -2
    a_33 = 3 - x_4
    a_34 = -x_3

    a_41 = 2*x_1
    a_42 = 2*x_2
    a_43 = 2*x_3
    a_44 = 0

    return [a_11 a_12 a_13 a_14; a_21 a_22 a_23 a_24; a_31 a_32 a_33 a_34; a_41 a_42 a_43 a_44]
end

function metodoGaussJordan(sistema, vectorResul)
    auxD = 0.0
    auxEscal = 0.0

    for i in 1:convert(Int64, size(sistema)[1])

        #Convirtiendo el valor x_ii en uno para escalonar
        auxD = sistema[i, i]
        for dividir in 1:convert(Int64, size(sistema)[1])
            sistema[i, dividir] = sistema[i, dividir] ./ auxD
        end
        vectorResul[i] = vectorResul[i] ./ auxD

        #Proceso para escalonar la matriz
        for j in i+1:convert(Int64, size(sistema)[1])
            auxEscal = sistema[j, i] #Obteniendo los valores

            for k in 1:convert(Int64, size(sistema)[1]) #Convirtiendo en cero los valores debajo de la diagonal principal
                sistema[j, k] = sistema[j, k] - (auxEscal * sistema[i, k])
            end
            vectorResul[j] = vectorResul[j] - (auxEscal * vectorResul[i])
        end
    end

    #Proceso inverso para convertir el sistema dado en una matriz diagonal
    for i in convert(Int64, size(sistema)[1]):-1:1

        for j in i-1:-1:1
            pivote = sistema[j, i] #Obteniendo los valores

            for k in convert(Int64, size(sistema)[1]):-1:1 #Convirtiendo en cero los valores arriba de la diagonal principal
                sistema[j, k] = sistema[j, k] - (pivote * sistema[i, k])
            end
            vectorResul[j] = vectorResul[j] - (pivote * vectorResul[i])
        end
    end

    #Regresa el valor del vector resultado obtenido de la resolución del sistema
    return vectorResul
end

function metodo_Newton(x_1n, x_2n, x_3n, x_4n)
    vector_xm = 0.0
    iteration = 0

    while iteration ≤ 1
        vectorFx = F(x_1n, x_2n, x_3n, x_4n)  #Vector resultante al evaluar X_n en la función
        matrizJ = J(x_1n, x_2n, x_3n, x_4n)   #Matriz Jacobiana evaluada en el vector X_n

        #Se obtiene el resultado equivalente de multiplicar
        #la matriz Jacobiana invertida multiplicada por el 
        #vector resultante al evaluar F(X_n), recordando que 
        #esto es lo mismo que resolver un sistema de ecuaciones 
        #J(X_n)*Y = F(X_n), en donde se busca al vector Y:
        #vectorResul = metodoGaussJordan(matrizJ, vectorFx)
        #Sin embargo, obtenemos el mismo resultado haciendo esto:
        vectorResul = inv(matrizJ)*vectorFx

        vector_xn = [x_1n, x_2n, x_3n, x_4n]

        vector_xm = vector_xn .- vectorResul    #Se obtiene a X_1 = X_0 - Y

        #Se iguala X_0 a X_1 para la siguiente iteración
        x_1n = vector_xm[1]
        x_2n = vector_xm[2]
        x_3n = vector_xm[3]
        x_4n = vector_xm[4]

        iteration = iteration + 1
    end

    return vector_xm
end

function sherman_morrison(x_1n, x_2n, x_3n, x_4n, x_1m, x_2m, x_3m, x_4m, A_0)
    x0 = [x_1n, x_2n, x_3n, x_4n]
    x1 = [x_1m, x_2m, x_3m, x_4m]
    y_i = F(x_1m, x_2m, x_3m, x_4m) .- F(x_1n, x_2n, x_3n, x_4n)
    s_i = x1 .- x0

    #Producto externo:
    AxY = A_0 * y_i
    prodExtern = 0.0
    for j in 1:convert(Int64, size(s_i)[1])
        prodExtern = prodExtern + s_i[j] * AxY[j]
    end

    A_1 = A_0 + ((s_i - (A_0 * y_i)) * (transpose(s_i) * A_0)) / prodExtern
    return A_1
end

function metodo_Broyden(x_1n, x_2n, x_3n, x_4n)
    iteration = 0
    δ = 1
    x_n = [x_1n, x_2n, x_3n, x_4n]

    #Primer iteración con Newton -> x_1
    x_m = metodo_Newton(x_1n, x_2n, x_3n, x_4n)

    mJacobi = inv(J(x_1n, x_2n, x_3n, x_4n)) #-> J(x_0)^-1

    iteration = 1

    while iteration ≤ max_iter && δ ≥ TOL

        mJacobi = sherman_morrison(x_n[1], x_n[2], x_n[3], x_n[4], x_m[1], x_m[2], x_m[3], x_m[4], mJacobi) #-> ≈ J(x_iteration)^-1
        x_n = [x_m[1], x_m[2], x_m[3], x_m[4]] #-> x_n = x_iteration
        
        #Se consigue a x_iteration+1
        x_m = x_n - (mJacobi * F(x_n[1], x_n[2], x_n[3], x_n[4]))

        #Error
        ω = abs.(x_m .- x_n)
        δ = max(ω[1], ω[2], ω[3], ω[4])

        iteration = iteration + 1
    end

    println("La raíz encontrada es X_$iteration = $x_m")
    println("Evaluando la función original con la raíz encontrada: ")
    resultadoF = F(x_m[1], x_m[2], x_m[3], x_m[4])
    println("F(X_$iteration) = $resultadoF")
end

#Valores iniciales

#Valores que acercan las dos primeras soluciones{
#=     x_1n = 1.0
    x_2n = 1.0
    x_3n = 1.0
    x_4n = 1.0  =#

#=     x_1n = -1.0
    x_2n = -1.0
    x_3n = -1.0
    x_4n = 1.0 =#
#}

#Valores que acercan las dos segundas soluciones{
#=     x_1n = -1
    x_2n = -0.5
    x_3n = 0.5
    x_4n = 3.5 =#

#=     x_1n = 1
    x_2n = 0.5
    x_3n = -0.5
    x_4n = 3.5 =#
#}

#Valores que acercan las dos segundas soluciones{
    x_1n = 1
    x_2n = -1
    x_3n = 1
    x_4n = 7

#=     x_1n = -1
    x_2n = 1
    x_3n = -1
    x_4n = 7 =#
#}

metodo_Broyden(x_1n, x_2n, x_3n, x_4n)