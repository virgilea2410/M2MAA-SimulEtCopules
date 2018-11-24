#Projet Simulation et Copules

#Vous trouverez dans un 1er temps le code pour l'affichage des courbes en 3 dimensions.
#Et ensuite l ensemble du script R de notre projet.

h = function
(x, y)
{
  return(1 + log(x + y) * 2 * x)
}

#Graphe de h
persp3d(h,
        xlim=c(0, 2),
        ylim=c(0, 3),
        main="fonction h",
        ylab="y",
        xlab="x")


# pour la fonction persp3d
require("rgl")
# pour la fonction quad2d
require("pracma")

#install.packages("rgl")

# nombres totales de simulations
n = 1000
# bornes de l'espace de définition de X
inf_x = 0
sup_x = 2
# bornes de l'espace de définition de Y
inf_y = 0
sup_y = 3
# En regardant la fonction dans persp3d, on peut voir que h est à valeur dans [0, 8]
inf_v = 0
sup_v = 8

H = quad2d(h, inf_x, sup_x, inf_y, sup_y)
paste("Vraie valeur de l'intégrale de h : ", H)


ux = runif(n, inf_x, sup_x)
uy = runif(n, inf_y, sup_y)

v = runif(n, inf_v, sup_v)

# Méthode A
# Marche bien mais on veut pouvoir afficher la convergence
#Sn = 1/n * sum(v <= h(ux, uy))
#I = Sn * ((sup_x-inf_x)*(sup_y-inf_y)*(sup_v-inf_v))
#paste("valeur de l'intégrale estimé par la méthode 1 : ", I)

# Méthode B
keeped_values = v <= h(ux, uy) * 1
Sn_array1 = 1/(1:n) * cumsum(keeped_values) * ((sup_x-inf_x)*(sup_y-inf_y)*(sup_v-inf_v))
Sn1 = 1/n * sum(keeped_values) * ((sup_x-inf_x)*(sup_y-inf_y)*(sup_v-inf_v))
biais1 = mean(Sn_array1[n]) - H
variance1 = var(Sn_array1 - H)



x = runif(n, inf_x, sup_x)
y = runif(n, inf_y, sup_y)

Sn_array2 = (sup_x-inf_x) * (sup_y-inf_y) * cumsum(h(x, y))/(1:n)
Sn2 =  1/n * (sup_x-inf_x) * (sup_y-inf_y) * sum(h(x, y))
biais2 = mean(Sn_array2[n]) - H
variance2 = var(Sn_array2 - H)


# On approxime l'intégral de f(x) par l'esperance de f(u) avec u une uniforme
# sur le support de x (si x défini sur [0, 3] alors u uniforme sur [0, 3])
# On applique la méthode des variables antithétiques : à chaque uniforme u, 
# on calcule l'esperance de 1/2 * (f(u) + f(sup(X) - u))
# pour l'exemple précedent, celà donne : 1/2 * (f(u) + f(3 - u))
# -> l'avantage principale est d'améliorer la vitesse de convergence :
# -> on a deux nouvelles observations (u et 1-u) à chaque nouvelle simulation de u 
# -> c'est la méthode des variables antithétiques
x = runif(n, inf_x, sup_x)
y = runif(n, inf_y, sup_y)

Sn_array21 = (sup_x-inf_x) * (sup_y-inf_y) * cumsum(h(x, y) + h(sup_x - x, sup_y - y))/(2*(1:n))
Sn21 =  1/n * (sup_x-inf_x) * (sup_y-inf_y) * sum(1/2 * (h(x, y) + h(sup_x - x, sup_y - y)))
biais21 = mean(Sn_array21[n]) - H
variance21 = var(Sn_array21 - H)


fx = function(x)
{
  return((0.2+0.3*x)*(x>=0 && x<=2))
}

fy = function(y)
{
  return((2/15+ 2/15*y)*(x>=0 && x<=3))
}

# on vérifie que les intégrales sur leurs supports respectifs
# valent 1
quad(fx, 0, 2)
quad(fy, 0, 3)


Fx = function(x)
{
  return((0.2*x+(0.15*x^2)))
}

Fy = function(y)
{
  return(2/15*x+ (1/30*x^2))
}


Fx_inv = function(x)
{
  return(sqrt(20/3 * x + 4/9) - 2/3)
}

Fy_inv = function(y)
{
  return(sqrt(15*y + 1) - 1)
}


#u = runif(n, inf_v, sup_v)
u = runif(n, 0, 1)
X = Fx_inv(u)
Y = Fy_inv(u)
Sn_array3 = 1/(1:n) * cumsum(h(X, Y)/(fx(X)*fy(Y)))
Sn3 = 1/n * sum(h(X, Y)/(fx(X)*fy(Y)))
biais3 = mean(Sn_array3[n]) - H
variance3 = var(Sn_array3 - H)


franck=function(u,v) {
  teta=0.55
  v=1-v
  return( (-teta*exp(-teta*u)*exp(-teta*v)*(exp(-teta)-1))/(((exp(-teta)-1)+(exp(-teta*u)-1)*(exp(-teta*v)-1))^2))}


h_copule=function(x,y) { fx(x)*fy(y)*franck(Fx(x),Fy(y))}

u = runif(n, 0, 1)
v=  runif(n, 0, 1)
X = Fx_inv(u)
Y = Fy_inv(v)

Sn_array4 = 1/(1:n) * cumsum(h(X,Y)/(h_copule(X,Y)))
Sn4 = 1/n * sum(h(X,Y)/h_copule(X,Y))
biais4 = mean(Sn_array3[n]) - H
variance4 = var(Sn_array3 - H)
# "Vrai valeur de l'intégrale de h
print(paste("Vraie valeur de l'intégrale de h : ", H),quote=FALSE)
# On affiche les résultat des approximations
print(paste("Méthode 1:"),quote=FALSE)
print(paste("Valeur de l'intégrale de h : ", Sn1,"biais: ",biais1, " variance", variance1 ),quote=FALSE) 
print(paste("Méthode 2:"),quote=FALSE)
print(paste("valeur de l'intégrale de h estimé par la méthode 2 : ", Sn2,"biais: ",biais2, " variance", variance2),quote=FALSE )
print(paste("Méthode 2.1:"),quote=FALSE)
print(paste("valeur de l'intégrale de h estimé par la méthode 2.1 : ", Sn21,"biais: ", biais21, " variance", variance21),quote=FALSE )
print(paste("Méthode 3:"),quote=FALSE)
print(paste("valeur de l'intégrale de h estimé par la méthode 3 : ", Sn3,"biais: ", biais3, " variance", variance3),quote=FALSE )
print(paste("Méthode 4:"),quote=FALSE)
print(paste("valeur de l'intégrale de h estimé par la méthode 4 : ", Sn4,"biais: ", biais4, " variance", variance4),quote=FALSE )


# On affiche le graph des 5 estimateurs afin d'étudier
# la vitesse de convergence vers la "vrai valeur"
plot(Sn_array1,
     col="gray0",
     type="l",
     ylim=c(10,20),
     #xlim=x(100, n),
     main="Vitesse de convergence",
     xlab="Nombres de simulations",
     ylab="Intégrale estimée")
abline(H, 0)
lines(Sn_array2,type="l", col="blue")
lines(Sn_array21, col="yellow")
lines(Sn_array3, col="green")
lines(Sn_array4, col="purple")




legend(400, 12.8,
       legend=c("Méthode des Volumes",
                "Méthode de Monte Carlo",
                "Méthode des Variables Antithétiques",
                "Méthode de la densité jointe",
                "Méthode des copules"),
       col=c("gray0",
             "blue",
             "yellow",
             "green",
             "purple"),
       lwd=4,
       cex=0.8)