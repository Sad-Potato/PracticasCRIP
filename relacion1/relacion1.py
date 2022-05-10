"""
    CRIPTOGRAFIA Y COMPUTACIÓN
    Relación 1: Aritmética modular

    @author: Sergio Fernández Vela

"""

from operator import mod
import numpy as np
import os
import random
import math
import time, collections

# Para comprobaciones
from sympy import legendre_symbol

######### Ejercicio 1 #########

def mcd_r(m):
    if m[0][0]%m[1][0]!=0:
        aux=(m[0][0]//m[1][0])
        c=[m[0][1]-m[1][1]*aux,m[0][2]-m[1][2]*aux]
        m[0][1:]=m[1][1:]
        m[1][1:]=c
        aux=m[1][0]
        m[1][0]=m[0][0]%m[1][0]
        m[0][0]=aux
        m=mcd_r(m)
    return m

def mcd(a,b):
    m=[[a,1,0],[b,0,1]]
    m=mcd_r(m)
    return m[1][:]


######### Ejercicio 2 #########

def modInverso(a,b):
    # Primero comprobamos que a y b sean primos relativos
    # es decir que el mcd de a y b sea 1

    # Calculamos con la función del apartado anterior
    # el mcd de a y b y sus coeficientes de bezout
    m=mcd(a,b)

    # Si el mcd es 1 continuamos si no, devolvemos
    # que no tiene inverso
    if m[0]==1:
        return m[1]%b
    else:
        return "Error: a no tiene inverso"

######### Ejercicio 3 #########

# Función para el calculo de a^b modulo n
def modPotencia(a,b,n):
    p=1
    while b>0:
        r=b%2
        if r==1:
            p=p*a%n
        a=(a*a)%n
        b=(b-r)//2
    return p

######### Ejercicio 4 #########

# Función para calcular los valores de 
# u y s para un número n-1 que es par,
# esto consiste en calcular la descomposición en 
# números primos de un número par

def descomposicionUyS(n):
    u=0
    while(modPotencia(n,1,2)==0):
        n=n//2
        u+=1
    return u,n

# Funciónes para hacer las comprobaciones pertinentes a 
# cada numero aleatorio que escogemos de Zₚ y sus potencias 
# de 2ᵘ*s, esta función se llama k veces desde "esPrimo" con 
# k dependiendo del grado de confianza que necesitemos

def listaL_apuntes(a,u,s,p):
    # Para cada a_i comprobamos unas condiciones
    # para sus potencias

    # Partimos de a^s
    aS=modPotencia(a,s,p)
    if aS==1 or aS==p-1:
        return True

    # Usamos u-1 ya que ya hemos comprobado 
    # para u=0
    for k in np.arange(u-1):
        aS=modPotencia(aS,2,p)
        # Si hay algun "a" que valga -1 devolvemos primo
        if aS==p-1:
            return True
        # Si hay algun "a" que valga 1 devolvemos que no es primo
        if aS==1:
            return False
    return False


# Función para determinar si dado un número p es (probablemente)
# primo usando el método de Miller-Rabin, si devolvemos que no es primo
# lo hacemos con total seguridad y si decimos que lo es lo hacemos 
# con un error de 1/4¹⁰ 

def esPrimo(p, k):

    # Primer caso, si p es par y p ≠ 2
    if p!=2 and p%2==0:
        return False

    # Escogemos k números aleatorios que 
    # estén en Zₚ y realizamos los siguientes
    # calculos para cada uno de ellos
    u,s=descomposicionUyS(p-1)

    for i in np.arange(k):
        # Escogemos el número aleatorio
        a_i=random.randint(2,p-1)

        # Llamamos a la función que se 
        # encarga de hacer las comprobaciones
        # pertinentes para el número en Zₚ
        # Si para alguno de los números aleatorios
        # devuelve que no es primo paramos y terminamos
        # ya que tenemos la seguridad de que no lo es
        if not listaL_apuntes(a_i,u,s,p):
            return False

    # Si aguantan las condiciones devolvemos
    # que es primo
    return True


######### Ejercicio 5 #########

# Para el cálculo de un s que sea mayor o igual que 
# √p usamos el método de Newton-Raphson para el calculo
# de raices cuadradas aunque sabemos la función 
# math.sqrt de python tambien funciona para valores grandes

# Valor n para el que queremos la raiz cuadrada y l
# para el nivel de precisión
def newtonRaphson(n, l):
    # Primero hacemos una suposición con
    # respecto al valor de la raiz cuadrada de n
    raiz=n

    # Condicion, iteramos mientras el error 
    # entre iteraciones sea mayor a l
    cond=True

    while(cond):

        raiz1=0.5*(raiz+(n/raiz))

        # Si hemos llegado a una raiz con un error menor a 
        # l paramos
        if(abs(raiz1-raiz)<l):
            cond=False
        
        # Actualizamos el valor de la raiz
        # con el calculado a partir de este
        raiz=raiz1

    return raiz


# Funcion para el algoritmo de paso enano-paso gigante
# El algoritm tiene como entrada un a y c no nulos, un p 
# primo y buscamos un b tal que a^b = c en Zₚ

def pasoEnanoGigante(a,c,p):
    
    # El primer paso es encontrar un "s" que 
    # sea mayor o igual que √p, lo calculamos 
    # con una precisión suficiente
    raizP=newtonRaphson(p,0.00000001)

    # El valor de s sera el siguiente entero
    # que exista a partir de la raiz de p
    # (s es un entero mayor o igual que raiz de p)
    s=math.ceil(raizP)
    
    # Calculamos la primera de las listas
    listaS={}
    elemento=c
    
    for r in np.arange(s): # Es decir [0,s-1]
        # Añadimos el valor y el indice a la tabla hash
        listaS[elemento]=r

        # Calculamos el valor de la 
        # siguiente iteración
        elemento=modPotencia(elemento*a,1,p)
        

    # A partir de la primera lista vamos calculando 
    # valores de la segunda hasta que encontremos 
    # una coincidencia o lleguemos al final

    # No necesitamos guardar los valores de la segunda lista
    # ya que estamos comprobandolos 1 a 1

    encontrado=False
    elemento=inicial=modPotencia(a,s,p)
    for t in np.arange(s)+1:

        # Paramos si el elemento calculado en la iteración actual
        # coincide con alguno de la primera lista
        if elemento in listaS:
            encontrado=True
            logaritmo=t*s-listaS[elemento] # r
            break

        # Calculamos el valor de la 
        # siguiente iteración
        elemento=modPotencia(elemento*inicial,1,p)

    # Si se ha encontrado alguna coincidencia entre las 2 listas
    # se devuelve el algoritmo calculado, en otro caso se vuelve
    # que no existe
    if encontrado:
        return logaritmo
    
    return "no existe el logaritmo"



######### Ejercicio 6 #########
"""
    Escribe una función que, dado un entero a y un primo p 
    con a p = 1, devuelve r tal que r 2 ≡ a mód p
    ([3, §2.3.4]; primero te hará falta implementar 
    el símbolo de Jacobi [1, 2.149]).
"""

# Función para el cálculo del simbolo de Jacobi para 
# un a sobre n (y para el cálculo del simbolo de Legendre en el 
# caso de que n sea primo) (a/p(n))

def simboloJacobiLegendre(a,n):
    # Comprobamos que se cumplan las condiciones
    # iniciales establecidas, aunque se puede 
    # calcular de igual manera aunque no se cumplan
    # estan condiciones
    assert n>=3 and a>=0 and a<n

    if a==0:
        return 0
    if a==1:
        return 1
    
    # Descomponemos a en potencias de 2 por un 
    # número impar
    e,a_1=descomposicionUyS(a)

    # Si el número es par entonces s=1
    s=0
    if modPotencia(e,1,2)==0:
        s=1
    else: 
        # Si no es impar comprobamos que n sea congruente con 1
        # modulo 8 o que lo sea con 7 en cuyo caso s=1
        if modPotencia(n-1,1,8)==0 or modPotencia(n-7,1,8)==0:
            s=1
        # Si no se cumplen las congruencias anteriores comprobamos que
        # n sea congruente con 3 o 5 modulo 8
        if modPotencia(n-3,1,8)==0 or modPotencia(n-5,1,8)==0:
            s=-1
    
    # Si n y a_1 son congruentes con 3 módulo 4 entonces 
    # s=-s
    if modPotencia(n-3,1,4)==0 and modPotencia(a_1-3,1,4)==0:
        s=-s
    
    n1=modPotencia(n,1,a_1)

    if a_1==1:
        return s
    else:
        return s*simboloJacobiLegendre(n1,a_1)
    

# A partir de la función para el cálculo del
# simbolo de Jacobi calculada procedemos con la 
# primera de las funciones pedidas

# La función tiene como entrada un entero a y
# un primo p 
def squaringTrapdoorRabin(a, p):

    # Ejemplo de output para el simbolo de Jacobi para la 
    # implementación propia y la implementación de sympy

    """ print([legendre_symbol(i, 7) for i in range(7)])
    print([simboloJacobiLegendre(i, 7) for i in range(7)]) """

    # A partir de la implementación anterior del 
    # símbolo de Jacobi tenemos como condición que el
    # símbolo de Legendre(ya que si p es primo el simbolo de
    # Jacobi es equivalente al de Legendre) para
    # a sobre p sea 1, ya que sabemos
    # que para que exista algún r tal que r^2 sea 
    # congruente con a modulo p el símbolo de Legendre 
    # para a y para p tiene que ser 1

    if simboloJacobiLegendre(a,p)!=1:
        print("No existe un r para el a y p dados")
        # Con devolver -1 nos basta, ya que en todos los demas
        # casos se devuelve un entero módulo p y me es util para
        # el apartado siguiente
        return -1
    
    # Si existe un r lo buscamos según el algoritmo
    # planteado por "Lecture Notes on Cryptography"
      
    # En primer lugar si p es congruente con 3 módulo 4, 
    # devolvemos que a elevado a m+1 módulo p es la raiz, 
    # siendo m el número por el que tenemos que multiplicar
    # 4 para que sea igual que p-3
    if modPotencia(p-3,1,4)==0:
        m=(p-3)//4
        return modPotencia(a,m+1,p)


    if modPotencia(p-1,1,4)==0:
        m=(p-1)//4
        
        # Buscamos un valor aleatorio de b que
        # satisfaga que su simbolo de Jacobi sobre
        # p sea -1, es decir que no sea un residuo
        # cuadrático
        b=random.randint(0,p-1)
        while simboloJacobiLegendre(b,p)!=-1:
            b=random.randint(0,p-1)

        # Cuando hemos encontrado un b
        i=modPotencia(2*m,1,p)
        j=0

        # Repetimos el siguiente bucle mientras que 
        # i sea par
        while modPotencia(i,1,2)!=1:
            i=i//2
            j=j//2

            a_i=modPotencia(a,i,p)
            b_i=modPotencia(b,j,p)
            if modPotencia(a_i*b_i,1,p)==p-1:
                j=modPotencia(j+2*m,1,p)
        
        # Cuando i pase a ser impar devolvemos que 
        # a elevado a i+1 entre 2 por b elevado a j
        # entre 2 es una raiz de a módulo p
        return modPotencia(modPotencia(a,(i+1)/2,p)*modPotencia(b,j/2,p),1,p)

# Función alternativa a la planteada arriba sacada de 
# "Handbook of Applied Cryptography" para sacar las
# 2 raices de un a módulo p

def raicesEnP(a, p):
    if not esPrimo(p,10):
        return "Error, "+str(p)+" no es primo "
    
    if modPotencia(p-3,1,4)==0:
        r=modPotencia(a,(p+1)/4,p)
        return [r,-r]
    
    if modPotencia(p-5,1,8)==0:
        d=modPotencia(a,(p-1)/4,p)

        if d==1:
            r=modPotencia(a,(p+3)/8,p)
            return [r,-r]
        if d==p-1:
            r=modPotencia(modPotencia(2*a,1,p)*modPotencia(4*a,(p-5)/8,p),1,p)
            return [r,-r]

# Segundo apartado del ejercicio, a partir de un "a" que es residuo
# cuadrático para un p y un q primos, uso el teorema chino de los 
# restos para calcular todas las raíces cuadradas de a modulo n siendo 
# n=pq a partir de las de p y q

# Uso el algoritmo 3.44 de "Handbook of applied Cryptography"

def raicesNCompuesto(a,p,q):
    n=p*q
    # Calculo las raices de "a" módulo p y q por separado
    rp,sq=squaringTrapdoorRabin(a,p),squaringTrapdoorRabin(a,q)
    rp2,sq2=modPotencia(-rp,1,p),modPotencia(-sq,1,q)

    if rp==-1 or sq==-1:
        return "Error, "+str(a)+" no es residuo cuadratico para "+str(p)+" o "+str(q)

    # Usamos el algoritmo extendido de euclides para encontrar c y d
    # tales que c*p+d*q=1 es decir los coeficientes de bezout de p y q
    m=mcd(p,q)
    c,d=m[1:]
    c,d=int(c),int(d)

    # Ahora con los coeficientes de Bezout calculamos las raices 
    # de n a partir de las de p y q usando el teorema chino de 
    # los restos
    x=modPotencia(rp*d*q+sq*c*p,1,n)
    y=modPotencia(rp*d*q+sq2*c*p,1,n)
    x2,y2=modPotencia(-x,1,n),modPotencia(-y,1,n)

    return [x,x2,y,y2]
    
    
######### Ejercicio 7 #########

# Todo número natural mayor que 1, o es primo, o se puede expresar
# de forma única como producto de primos

# Factorización de un número por el método de 
# Fermat, metodo útil en el caso de que n tenga dos 
# divisores relativamente proximos y proximos a 
# n/2
def factorizacionFermat(n):

    # Obtenemos la raiz cuadrada de n y obtenemos el entero
    # inmediatamente superior, mientras x^2-n no sea un 
    # cuadrado perfecto incrementamos x en una unidad

    x=math.ceil(newtonRaphson(n,0.000001))
    x2n=pow(x,2)-n

    # pow(round(newtonRaphson((pow(x,2)-n),0.00001)),2)!=(pow(x,2)-n):
    while pow(math.floor(newtonRaphson(x2n,0.0001)+0.5),2)!=x2n:
        # Compruebo que sea un cuadrado perfecto calculando la aproximacion
        # de la raiz cuadrada y en el caso de que el cuadrado de esta apro-
        # mación sea igual pues concluyo que lo es
        x+=1
        x2n=pow(x,2)-n
    
    # Una vez que tenemos un x tal que x^2-n es un cuadrado perfecto
    # le asignamos a y el valor de la raiz cuadrada de x^2-n
    y=int(newtonRaphson((pow(x,2)-n),0.00001)+0.5)

    # Devolvemos los factores encontrados
    return [x+y,x-y]
    
# Factorización de un número por el método de 
# Pollard, devolvemos un factor d no trivial de n
def factorizacionPollard(n):
    # Comprobamos que sea primo, en cuyo caso devolvemos 
    # el propio n
    if esPrimo(n,10):
        return n
    
    a,b,d=2,2,1
    # Mientras d valga 1
    while d==1:
        
        # Le aplicamos a "a" la funcion g(x)=(x^2+1) mod n 1 vez
        # y a b dos veces 
        a=modPotencia(modPotencia(a,2,n)+1,1,n)
        b=modPotencia(modPotencia(b,2,n)+1,1,n)
        b=modPotencia(modPotencia(b,2,n)+1,1,n)

        # Calculamos el MCD para a-b y n y se
        # lo asignamos a d
        d=mcd(modPotencia(a-b,1,n),n)[0]
    
        if 1<d<n:
            return d
        if d==n:
            return "Fracaso al encontrar un factor no trivial"
    


######### Ejercicio 8 #########

# Benchmark de tiempos acordada entre
# compañeros de la asignatura
def tiempos(primos):
    tiempos = collections.defaultdict(int)
    EJECUCIONES = 1000
    
    for _ in range(EJECUCIONES):
        t = time.time()
        mcd(primos[0], primos[1])
        tiempos['gcd'] += time.time() - t
        
        t = time.time()
        modInverso(primos[0], primos[1])
        tiempos['modinv'] += time.time() - t
        
        t = time.time()
        modPotencia(primos[0], primos[0], primos[1])
        tiempos['modpow'] += time.time() - t
        
        t = time.time()
        esPrimo(primos[0],10)
        tiempos['Miller-Rabin'] += time.time() - t
        
        aux = modPotencia(51, 79, primos[0])
        t = time.time()
        pasoEnanoGigante(51, aux, primos[0])
        tiempos['enano-gigante'] += time.time() - t
        
        aux = modPotencia(primos[0], 2, primos[1])
        t = time.time()
        squaringTrapdoorRabin(aux, primos[1])
        tiempos['Tonelli-Shanks'] += time.time() - t
        
        n = primos[0] * primos[1]
        aux = modPotencia(57, 2, n)
        t = time.time()
        raicesNCompuesto(aux, primos[0], primos[1])
        tiempos['CRT'] += time.time() - t
        
        t = time.time()
        factorizacionPollard(n)
        tiempos['Pollard'] += time.time() - t
 
    for i in tiempos: tiempos[i] /= EJECUCIONES
 
    # Fermat va mucho más lento, así que lo hago con menos ejecuciones
    for _ in range(25):
        t = time.time()
        factorizacionFermat(n)
        tiempos['Fermat'] += time.time() - t
    tiempos['Fermat'] /= 25


    for i in tiempos:
        print(i+':', tiempos[i]*1000000)

"""
    Los tiempos que he obtenido son los siguientes, en microsegundos: 

                Tiempos personales          | Comparativa de un compañero
    Ejercicio 1     8.23µs                      7.316µs
    Ejercicio 2     8.32µs                      7.515µs
    Ejercicio 3     2.77µs                      7.832µs
    Ejercicio 4     33.49µs                     93.57µs
    Ejercicio 5     67.89µs                     543µs
    Ejercicio 6.1   15.68µs                     171µs
    Ejercicio 6.2   61.96µs                     317µs
    Ejercicio 7.1   558027µs(0.558 seg)         1225000µs             
    Ejercicio 7.2   2533.74µs(2.533 ms)         1937µs                                  


"""


###########################################################
# Main para elegir el ejercicio que queremos mostrar
###########################################################

def main():
    clear = lambda: os.system('clear')
    tecla=""
    random.seed(42)

    #################################################################################
    # PARAMETROS PARA LOS EJERCICIOS
    #################################################################################

    a,b=28,13
    c,d=252336560693540533935881068298825202079,38942750026936412998460304986028600003
    e,f,g=13,2,19
    # Lista de primos de gran tamaño para hacer 
    # pruebas
    primos=[46381, 768479, 9476407, 36780481, 562390847, 1894083629,
    65398261921, 364879542899, 8590365927553, 28564333765949, 123456789101119,
    623084000430982607975364879776457600049]

    while(tecla!="e"):
        print("Elija un ejercicio (1-8)")
        tecla=input()
        #clear() # Limpiamos la terminal

        #################################################################################

        if tecla=='1':
            m=mcd(a,b)
            print("El MCD positivo de",a,"y "+str(b)+" = "+str(m[0])+" ;u:",m[1],"v:",m[2])
        if tecla=='2':
            print("El inverso de",a,"es",modInverso(a,b))
            """ print(modPotencia(7*28,1,b)) """
        if tecla=='3':
            # Probamos nuestra función con números muy grandes, comprobado en gap
            print("El número",c,"elevado a",d,"módulo",primos[-2],"es igual a:")
            print(modPotencia(c,d,primos[-2]))
        if tecla=='4':
            # Probamos miller rabin para primos muy grandes
            print("Lista de números primos:")
            for p in primos:
                res=": no es primo"
                if esPrimo(p,10):
                    res=": es primo"
                print("\t",p,res)
        
        if tecla=='5':
            print("El logaritmo en base",e,"de",f,"modulo",g,"es: ",
                pasoEnanoGigante(e,f,g))

        if tecla=='6':
            # Primer apartado
            a,p=6,53
            r=squaringTrapdoorRabin(a,p)
            print("Para",a,"y",p,"tenemos que",
                r,"es una raiz de",
                a,"módulo",p,"ya que",r,"^2 congruente con",a,"módulo",p,"es",
                modPotencia(modPotencia(r,2,p)-a,1,p))

            # Segundo apartado
            a,p,q=4,7,11
            n=p*q
            raices=raicesNCompuesto(a,q,p)
            print("\nLas raices de",n,"siendo n",p,"por",q,"son",raices)

            # Comprobacion de que para las 4 raices si se elevan al cuadrado
            # y se les resta a módulo n el resultado es 0
            """ print(modPotencia(modPotencia(raices[0],2,n)-a,1,n),
                modPotencia(modPotencia(raices[1],2,n)-a,1,n),
                modPotencia(modPotencia(raices[2],2,n)-a,1,n),
                modPotencia(modPotencia(raices[3],2,n)-a,1,n)) """

        if tecla=='7':
            n=6352351
            fact=factorizacionFermat(n)
            print("La factorización del número",n,"por Fermat es",fact[0],"y",fact[1])

            n2=455459
            fact2=factorizacionPollard(n2)
            print("Un factor no trivial obtenido de la factorización por Pollard de",n,"es",fact2)

        if tecla=='8':
            tiempos(primos)

        if tecla=='9':
            print(mcd(9,5))


        #################################################################################

        tecla=input("\nCualquier tecla para volver al menu, \"e\" para salir\n")

if __name__ == "__main__":
    main()
