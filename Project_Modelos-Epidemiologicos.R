
# Este programa hace uso de la librería "deSolve" para la resolución de las 
# respectivas ecuaciones diferenciales

######## Significado de los parámetros a recibir por las respectivas funciones ##########

# beta    = Tasa de infeccion
# gamma   = Tasa de recuperacion
# mu      = Tasa de muerte
# delta   = Tasa de nacimientos
# Estados = Vector con información sobre el número Personas o población que es suceptible a la infección,
#           personas infectadas, recuperadas etc.
#           debe ser llenado en orden SIR [Susceptibles, Infectados, Recuperados] o como corresponda
# sigma   = tasa de exposicion

########################## Modelo Epidemológico SIR ##############################
#install.packages("deSolve")
library(deSolve)
library(grid)
library(ggplot2)
library(dplyr)
SIR(Tiempo=50,Estados=c(1500,1,15),  beta=2.5,gamma=.1, delta = 1,mu=.01,Grafica = 1)
SEIR(Tiempo=20, Estados=c(100,15,25,1),beta=2.5,gamma=1, delta = 1,mu=1,Grafica = 1,sigma=.005)
SIS(Tiempo=20, Estados=c(100,5), beta=.01,gamma = .001, Grafica =1)

Zombieland(Tiempo=60, Estados=c(100,5,0), beta=.08, alfa=1, dseta=1, Grafica = 1)
  

SIR <- function(Tiempo, Estados, beta, gamma, delta = 0, mu = 0, Grafica = 0){
  if(length(Tiempo) == 1){                  # El parámetro "Tiempo" puede ser un vector o un sólo número, de ser un 
    Tiempo <- c(0:Tiempo)                   # dígito, se asignará un vector empezando en Tiempo[1] = 0
  }
  
  Condiciones <- c(Susceptibles = Estados[1],
                   Infectados   = Estados[2],
                   Recuperados  = Estados[3],
                   Poblacion    = sum(Estados)
  )

  Parametros <- c(B = beta, Y = gamma,
                  D = delta, M = mu
                )
                                            # Resolución de las ecuaciones diferenciales del modelo, función
                                            # llamada posteriormente con "ode" de la librería "deSolve"
  
  CalDif <- function(Tiempo, Estados, Parametros) {
    with(as.list(c(Estados, Parametros)), {
      dS <- (-1) * B * Susceptibles * Infectados + D - M * Susceptibles
      dI <-  B * Susceptibles * Infectados - Y * Infectados - M * Infectados
      dR <-  Y * Infectados - M * Recuperados
      Poblacion <- dS + dI + dR
      return(list(c(dS, dI, dR, Poblacion)))
    })
  }
  
  Datos <- ode(y = Condiciones,             # Obtención de resultados de las ecuaciones diferenciales ordinarias resueltas
             times = Tiempo,
             func = CalDif,
             parms = Parametros
  )   
  
  if(!delta){                               # De no ser agregada alguna tasa para dinámica vital, se descarta la
    Datos <- Datos[1:length(Datos[,1]),1:4] # información con respecto al número de población, esta se mantiene constante
  }

  if(Grafica){
    plot(Tiempo ,Datos[,2], lwd = 2,
         col = "orange", type = "l",
         main = "Modelo SIR",
         ylab = "Población",
         ylim = c(0, Datos[1,2] + Datos[1,3] + Datos[1,4])
    )
    grid(lty = "dotted")
    lines(Tiempo, Datos[,3], lwd = 2, col = "blue")
    lines(Tiempo, Datos[,4], lwd = 2, col = "black")
    legend(x = 0,y = -sum(Estados)/3, c("Susceptibles","Infectados", "Recuperados"),
           cex = 0.7, fill = c("orange", "blue", "black"), bty = "n", xpd = T
    )
  }
  return(Datos)
}

################## Modelo Epidemológico SEIR ########################


SEIR <- function(Tiempo, Estados, beta, gamma, sigma, mu= 0, delta = 0, Grafica = 0){
  
  if(length(Tiempo) == 1){
    Tiempo <- c(0:Tiempo)
  }
  
  Condiciones <- c(Susceptibles = Estados[1],
                   Exposicion   = Estados[2],
                   Infectados   = Estados[3],
                   Recuperados  = Estados[4],
                   Poblacion    = sum(Estados)
  )
  
  Parametros <- c(B = beta, Y = gamma, Q = sigma, M = mu, D = delta)
  
  CalDif <- function(Tiempo, Estados, Parametros) {
    with(as.list(c(Estados, Parametros)), {
      dS <-  D * Poblacion - B * Susceptibles * Infectados - M * Susceptibles
      dE <-  B * Susceptibles * Infectados - Q * Exposicion - M * Exposicion 
      dI <-  Q * Exposicion - Y * Infectados - M * Infectados
      dR <-  Y * Infectados - M * Recuperados
      Poblacion <- dS + dE + dI + dR
      return(list(c(dS ,dE, dI, dR, Poblacion)))
    })
  }
  
  Datos <- ode(y = Condiciones,
               times = Tiempo,
               func = CalDif,
               parms = Parametros
  )
  
  if(!delta){
    Datos <- Datos[1:length(Datos[,1]),1:5]
  }
  
  if(Grafica){
    plot.new()
    plot(Tiempo ,Datos[,2], lwd = 2,
         col = "orange", type = "l",
         main = "Modelo SEIR",
         ylab = "Población",
         ylim = c(0 , Datos[1,2] + Datos[1,3] + Datos[1,4] + Datos[1,5])
    )
    grid(lty = "dotted")
    lines(Tiempo, Datos[,3], lwd = 2, col = "blue")
    lines(Tiempo, Datos[,4], lwd = 2, col = "black")
    lines(Tiempo, Datos[,5], lwd = 2, col = "red")
    legend(x = 0, y = -sum(Estados)/3,
           c("Susceptibles", "Expuestos","Infectados", "Recuperados"),
           cex = 0.65,
           fill = c("orange", "blue", "black", "red"), 
           xpd = TRUE, bty = "n"
    )
  }
  return(Datos)
}

################# Modelo Epidemológico SIS #################### 

SIS <- function(Tiempo, Estados, beta, gamma, Grafica = 0){
  
  if(length(Tiempo) == 1){
    Tiempo <- c(0:Tiempo)
  }
  
  Condiciones <- c(Susceptibles = Estados[1] ,
                   Infectados   = Estados[2]
  )
  
  Parametros <- c(B = beta, Y = gamma)
  
  CalDif <- function(Tiempo, Estados, Parametros) {
    with(as.list(c(Estados, Parametros)), {
      dS <- (-1) * B * Susceptibles * Infectados + Y * Infectados
      dI <-  B * Susceptibles * Infectados - Y * Infectados
      return(list(c(dS, dI)))
    })
  }
  
  Datos <- ode(y = Condiciones,
               times = Tiempo,
               func = CalDif,
               parms = Parametros
  )
  
  if(Grafica){
    plot.new()
    plot(Tiempo ,Datos[,2], lwd = 2,
         col = "orange", type = "l",
         main = "Modelo SIS",
         ylab = "Población",
         ylim = c(0, Datos[1,2] + Datos[1,3])
    )
    grid(lty = "dotted")
    lines(Tiempo, Datos[,3], lwd = 2, col = "blue")
    legend(x = 0, y = -sum(Estados)/3, c("Susceptibles","Infectados"),
           cex = 0.9, fill = c("orange", "blue"),
           xpd = TRUE, bty = "n"
    )
  }
  return(Datos)
}

##################### ¡¡¡ ZOMBIES !!! #########################

####### Parámetros #######

# alfa  = Tasa de "Muerte" zombie
# dseta = Tasa en que humanos muertos se vuelven zombies
# beta  = Tasa de transmisión

Zombieland <- function(Tiempo, Estados, beta, alfa, dseta, Grafica = 0){
  
  if(length(Tiempo) == 1){
    Tiempo <- c(0:Tiempo)
  }
  
  Condiciones <- c(Susceptibles = Estados[1] ,
                   Zombies      = Estados[2] ,
                   Muertos      = Estados[3]
  )
  
  Parametros <- c(B = beta, A = alfa, Z = dseta)
  
  CalDif <- function(Tiempo, Estados, Parametros) {
    with(as.list(c(Estados, Parametros)), {
      dS <- (-1) * B * Susceptibles * Zombies
      dZ <-  (B - A) * Susceptibles * Zombies + Z * Muertos
      dM <- A * Susceptibles * Zombies - Z * Muertos
      return(list(c(dS, dZ, dM)))
    })
  }
  
  Datos <- ode(y = Condiciones,
               times = Tiempo,
               func = CalDif,
               parms = Parametros
  )
  
  if(Grafica){
    plot.new()
    plot(Tiempo ,Datos[,2], lwd = 2,
         col = "orange", type = "l",
         main = "¡¡¡ Zombies !!!",
         ylab = "Población",
         ylim = c(0, Datos[1,2] + Datos[1,3] + Datos[1,4])
    )
    grid(lty = "dotted")
    lines(Tiempo, Datos[,3], lwd = 2, col = "blue")
    lines(Tiempo, Datos[,4], lwd = 2, col = "black")
    legend(x = 0, y = -sum(Estados)/3, c("Susceptibles","Zombies", "Muertos"),
           cex = 0.8, fill = c("orange", "blue", "black"),
           xpd = TRUE, bty = "n"
           )
  }
  return(Datos)
}


##### Gif para ver evolución de infección #####

#### Este gif acepta números enteros

GIFFun <- function(infectados,poblacion){
  num_renglones <- sqrt(poblacion)
  png(file="Poexample%02d.png", width=200, height=200)
  for (ni in 1:length(infectados)){
    grid.newpage()
    contador<-1;
    x <- sample(1:poblacion,infectados[ni])
    matriz_device <- grid.layout(nrow =num_renglones + 1 , ncol = num_renglones)
    pushViewport(viewport(layout = matriz_device))
    for (i in 1: nrow(matriz_device)){                  
      for (j in 1:ncol(matriz_device)){                 
                                                          
        pushViewport(viewport(layout.pos.row= i, layout.pos.col= j))
        if(length(x[x==contador])){
          grid.circle(x=0.5,y=0.5,gp=gpar(fill="red",col="black"))
        }
        else{grid.circle(x=0.5,y=0.5,gp=gpar(fill="Green",col="black"))}
        popViewport()
        if(contador==poblacion){break}
        contador<-contador+1
      }
    }
    
    showViewport(viewport(layout = matriz_device))
    ##break
  }
  dev.off()
  system("convert -delay 80 *.png example_1.gif")
  file.remove(list.files(pattern=".png"))
}
