"""
    CRIPTOGRAFIA Y COMPUTACIÓN
    Relación 1: Aritmética modular

"""

from operator import mod
import numpy as np
import os
import random
import math

from sympy import legendre_symbol

######### Ejercicio 1 #########

def mcd_r(m):
    if m[0,0]%m[1,0]!=0:
        c=m[0,1:]-(m[0,0]//m[1,0])*m[1,1:]
        m[0,1:]=m[1,1:]
        m[1,1:]=c
        aux=m[1,0]
        m[1,0]=m[0,0]%m[1,0]
        m[0,0]=aux
        m=mcd_r(m)
    return m

def mcd(a,b):
    m=np.array([[a,1,0],[b,0,1]])
    m=mcd_r(m)
    return m[1,:]


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
    a_s=modPotencia(a,s,p)
    if a_s==1 or a_s==p-1:
        return True

    # Usamos u-1 ya que ya hemos comprobado 
    # para u=0
    for k in np.arange(u-1):
        a_s=modPotencia(a_s,2,p)
        if a_s==p-1:
            return True
        if a_s==1:
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
    # iniciales
    """ assert n>=3 and a>=0 and a<n """

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
        return "No existe un r para el a y p dados"
    
    # Si existe un r lo buscamos

    for m in np.arange(p-1)+1:
        
        if modPotencia(4*m+3,1,p)==0:
            print("#### 1 ####")
            print(a,m+1,p," #### ",modPotencia(-1,(p-1)/2,p))
            return modPotencia(a,m+1,p)
    
        if modPotencia(4*m+1,1,p)==0:
            
            # Buscamos un valor aleatorio de b que
            # satisfaga que su simbolo de Jacobi sobre
            # p sea -1
            b=random.randint(0,p-1)
            while simboloJacobiLegendre(b,p)==-1:
                b=random.randint(0,p-1)

            i=modPotencia(2*m,1,p)
            j=0

            while modPotencia(i,1,2)!=0:
                i=i/2
                j=j/2
                if modPotencia(modPotencia(a,i,p)*modPotencia(b,j,p),1,p)==-1:
                    j=modPotencia(j+2*m,1,p)
            
            print("#### 2 ####")
            return modPotencia(modPotencia(a,(i+1)/2,p)*modPotencia(b,j/2,p),1,p)
   

    """ for m in np.arange(p):
        
        if modPotencia(4*m+3,1,p)==0:
            print("#### 1 ####")
            return int(modPotencia(a,m+1,p))
    for m in np.arange(p):
        if modPotencia(8*m+5,1,p)==0 and modPotencia(modPotencia(a,2*m+1,p)-1,1,p)==0:
            print("#### 2 ####")
            return int(modPotencia(a,m+1,p))
    for m in np.arange(p):
        if modPotencia(8*m+5,1,p)==0 and modPotencia(modPotencia(a,2*m+1,p)+1,1,p)==0:
            print("#### 3 ####")
            return int(modPotencia(0.5*modPotencia(4*a,m+1,p)*(p+1),1,p))
     """



    """     u,s=descomposicionUyS(p-1)
    r=modPotencia(a,(s+1)/2,p)

    a_s=modPotencia(modPotencia(r,2,p)*modPotencia(a,p-1,p),1,p)

    m=random.randint(0,p-1)
    while simboloJacobiLegendre(m,p)!=-1:
        m=random.randint(0,p-1)
    
    mu=modPotencia(m,s,p) """
    
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
    while pow(round(newtonRaphson((pow(x,2)-n),0.00001)),2)!=(pow(x,2)-n):
        # Compruebo que sea un cuadrado perfecto calculando la aproximacion
        # de la raiz cuadrada y en el caso de que el cuadrado de esta apro-
        # mación sea igual pues concluyo que lo es
        x+=1
    
    # Una vez que tenemos un x tal que x^2-n es un cuadrado perfecto
    # le asignamos a y el valor de la raiz cuadrada de x^2-n
    y=round(newtonRaphson((pow(x,2)-n),0.00001))

    return [x+y,x-y]
    

def factorizacionPollard(n):


###########################################################
# Main para elegir el ejercicio que queremos mostrar
###########################################################

def main():
    clear = lambda: os.system('clear')
    tecla=""
    random.seed(42)

    while(tecla!="e"):
        print("Elija un ejercicio (1-8)")
        tecla=input()
        #clear() # Limpiamos la terminal

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

        #################################################################################

        if tecla=='1':
            m=mcd(a,b)
            print("MCD positivo de",a,"y "+str(b)+"="+str(m[0])+" ;u:",m[1],"v:",m[2])
        if tecla=='2':
            print(modInverso(a,b))
        if tecla=='3':
            # Probamos nuestra función con números muy grandes
            print(modPotencia(c,d,primos[-1]))
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
            a_s=[1,3,4,5,9,12,14,15,16,20,23]
            a_s=[4,6,7,9,10,11,13]
            p=53
            for a in a_s:
                r=squaringTrapdoorRabin(a,p)
                print(modPotencia(modPotencia(r,2,p)-a,1,p)," - ",modPotencia(r,2,p))
        
        if tecla=='7':
            n=6352351
            fact=factorizacionFermat(n)
            print("La factorización del número",n,"por Fermat es",fact[0],"y",fact[1])

            n2=7011461
            fact2=factorizacionPollard(n2)
            print("La factorización del número",n,"por Pollard es",fact2[0],"y",fact2[1])


        #################################################################################

        tecla=input("\nCualquier tecla para volver al menu, \"e\" para salir\n")

if __name__ == "__main__":
    main()
