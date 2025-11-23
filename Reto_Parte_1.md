Reto: Modelo SIR y Vacunación Parte 1
================
Aldo Resendiz
Noviembre de 2025

- [El modelo SIR](#el-modelo-sir)
- [Gráficos de la evolución del
  sistema](#gráficos-de-la-evolución-del-sistema)
- [Pregunta 1](#pregunta-1)
- [Pregunta 2 Reducción de
  Susceptibles](#pregunta-2-reducción-de-susceptibles)

## El modelo SIR

Consideremos un modelo para describir la dinámica de un grupo de
individuos de una población con exposición a una enfermedad que puede
contagiarse entre los miembros de la población. Esto puede modelarse
como un sistema dinámico denominado $`SIR`$ para una población de $`N`$
individuos en la que se considera la interacción entre un conjunto de
$`S`$ individuos *suceptibles* de contraer la enfermedad, un conjunto
$`I`$ de individuos *infectados* y uno conjunto $`R`$ de individuos
*recuperados* de la enfermedad.

Este modelo tiene los siguientes supuestos:

- la probabilidades de infectarse son iguales para todos los individuos
  de la población;

- la población es homogénea, es decir que los riesgos de infectarse son
  iguales para toos los suceptibles y que los tiempos para recuperarse
  son iguales para todos los infectados; y

- el tamaño $`N`$ de la población es constante.

El modelo maneja los diferentes conjuntos $`S`$, $`I`$ y $`R`$ como si
fueran compartimentos bien separados y considera que los individuos
pueden pasr de uno a otro en el caso de que se enfermen (cambio
$`S\rightarrow I`$) o que una vez enfermos se recuperen (cambio
$`I\rightarrow`$). Ademas, se asume que un individuo no puede pasar del
conjunto de suceptibles directamente al conjunto de recuperados.

Con estos supuestos y consideraciones, las ecuaciones diferenciales del
modelo SIR son:
``` math

\begin{aligned}
\frac{dS}{dt}&= -\beta \frac{I}{N} S\\
\frac{dI}{dt}&= \beta\frac{I}{N}S-\gamma I\\\
\frac{dR}{dt}&= \gamma I
\end{aligned}
```
donde:

- N=S+R+I

- la cantidad $`\beta\frac{I}{N}`$ representa la razón con que las
  personas salen del compartimento S (se infectan);

- en la primera ecuación $`dS`$ representa el cambio debido a las
  personas que salen del compartimento $`S`$ (el signo negativo se debe
  a que las personas salen)

- en la segunda ecuación $`dI`$ representa el cambio debido a las
  personas que salen del compartimento $`I`$ (una parte se debe a las
  personas que del compartimento $`S`$ pasan al compartimento $`I`$, y
  otra parte se debe a las personas que salen del compartimento $`I`$
  porque se recuperan);

- la cantidad $`\gamma`$ representa la razón con que las personas se
  recuperan.

``` r
# PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)

initial_state_values <- c(S = 999999, I = 1, R = 0)
parameters <- c(beta = 1, gamma = 0.1)
times <- seq(from = 0, to = 60, by = 1)   

sir_model <- function(time, state, parameters) {  
    with(as.list(c(state, parameters)), {
        N <- S+I+R 
        lambda <- beta * I/N
        dS <- -lambda * S                
        dI <- lambda * S - gamma * I   
        dR <- gamma * I                  
        return(list(c(dS, dI, dR))) 
    })
}

output <- as.data.frame(ode(y = initial_state_values, 
                            times = times, 
                            func = sir_model,
                            parms = parameters))
```

## Gráficos de la evolución del sistema

``` r
# Gráfico del modelo base
output_long <- melt(as.data.frame(output), id = "time")                  

ggplot(data = output_long, aes(x = time, y = value, colour = variable)) +  
  geom_line(linewidth = 1) +                                           
  xlab("Tiempo (días)")+                                         
  ylab("Número de individuos") +                                       
  labs(title = "Modelo SIR Básico", colour = "Subconjunto") +
  theme_minimal() +
  theme(legend.position = "bottom")
```

![](Reto_Parte_1_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

Con el modelo SIR se define la constante
``` math
R_0=\frac{\beta}{\gamma}
```
que representa el número de personas que cada contagiado infecta. Para
que la enfermedad analizada logre dispararse en forma de una epidemia
debe cumplirse que $`R_0 > 1`$.

También se define
``` math
R_{eff}=R_0\frac{S}{N}
```
que corresponde al número promedio de personas que cada contagiado
infecta. Este segundo valor $`R_{eff}`$ toma en cuenta de que durante la
evolución de la pandemia, al aumentar del número de personas inmunes en
la población cada persona contagiada infectará a un número de personas
cada vez menor.

## Pregunta 1

Analizando el dataframe “output” encuentre el día en que el número de
contagios es máximo.

``` r
# Encontrar la fila con el valor máximo de infectados (I)
fila_max_I <- output[which.max(output$I), ]
dia_pico <- fila_max_I$time
max_infectados <- fila_max_I$I

cat("El pico de contagios ocurre en el día:", dia_pico, "\n")
```

    ## El pico de contagios ocurre en el día: 18

``` r
cat("Número máximo de infectados:", round(max_infectados))
```

    ## Número máximo de infectados: 669741

*Relación analítica en el máximo*:Para encontrar la relación entre los
parámetros en el máximo de la curva de infección ($`I`$), analizamos la
ecuación diferencial de $`dI/dt`$:
``` math
 \frac{dI}{dt} = \beta \frac{I}{N} S - \gamma I 
```
En el punto máximo (el pico), la pendiente es cero ($`dI/dt = 0`$).
Entonces:
``` math
 0 = I \left( \beta \frac{S}{N} - \gamma \right) 
```
Como $`I \neq 0`$, debemos tener:
``` math
 \beta \frac{S}{N} = \gamma \implies S = \frac{\gamma}{\beta} N 
```
O equivalentemente, usando el número reproductivo básico
$`R_0 = \beta / \gamma`$:
``` math
 S_{pico} = \frac{N}{R_0} 
```
Esto significa que el pico ocurre exactamente cuando la población
susceptible $`S`$ cae hasta alcanzar el umbral $`N/R_0`$.

## Pregunta 2 Reducción de Susceptibles

Analizando el dataframe encuentre después de cuántos días el número de
“susceptibles” se reduce a la mitad.

``` r
S_inicial <- initial_state_values["S"]
umbral_mitad <- S_inicial / 2

# Filtrar los tiempos donde S ya es menor o igual a la mitad
fila_mitad <- subset(output, S <= umbral_mitad)[1, ] # Tomamos el primer dato

cat("El número de susceptibles cae a la mitad aprox. en el día:", fila_mitad$time, "\n")
```

    ## El número de susceptibles cae a la mitad aprox. en el día: 16

``` r
cat("Valor de S en ese día:", round(fila_mitad$S))
```

    ## Valor de S en ese día: 353135

Fórmula analítica (Aproximación SI): Si asumimos que al inicio la
epidemia crece muy rápido y ignoramos momentáneamente la recuperación
(comportamiento tipo logístico $`SI`$), la ecuación para $`S`$ es
aproximadamente:
``` math
 \frac{dS}{dt} \approx -\beta \frac{S(N-S)}{N} 
```
La solución analítica para $`S(t)`$ es la curva logística inversa:
``` math
 S(t) = \frac{N}{1 + \frac{I_0}{S_0} e^{\beta t}} 
```
Queremos hallar $`t`$ cuando $`S(t) = N/2`$. Sustituimos y despejamos:
``` math
 \frac{N}{2} = \frac{N}{1 + \frac{I_0}{S_0} e^{\beta t}} \implies 2 = 1 + \frac{I_0}{S_0} e^{\beta t} \implies 1 = \frac{I_0}{S_0} e^{\beta t} 
```
Tomando logaritmo natural:
``` math
 \ln(1) = \ln\left(\frac{I_0}{S_0}\right) + \beta t \implies 0 = \ln(I_0) - \ln(S_0) + \beta t 
```
``` math
 t = \frac{\ln(S_0) - \ln(I_0)}{\beta} 
```
\## Pregunta 3 Variación de $`\beta`$ (Fuerza de infección) Estudiamos
la dinámica manteniendo $`\gamma = 0.1`$ y variando $`\beta`$.Valores de
$`\beta`$: 0.1, 0.3, 0.7, 0.9, 1.2.

``` r
betas <- c(0.1, 0.3, 0.7, 0.9, 1.2)
gamma_fijo <- 0.1
resultados_p3 <- list()

for(b in betas){
  pars <- c(beta = b, gamma = gamma_fijo)
  out <- as.data.frame(ode(y=initial_state_values, times=times, func=sir_model, parms=pars))
  out$beta <- as.factor(b)
  resultados_p3[[length(resultados_p3)+1]] <- out
}

data_p3 <- do.call(rbind, resultados_p3)

ggplot(data_p3, aes(x=time, y=I, color=beta)) +
  geom_line(size=1) +
  labs(title = "Dinámica variando Beta (Gamma = 0.1)", 
       y = "Infectados", color = "Beta") +
  theme_minimal()
```

![](Reto_Parte_1_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->
*Relación necesaria*: Para que ocurra epidemia ($`dI/dt > 0`$ al
inicio), se necesita que $`\beta > \gamma`$.Esto define el Número
Reproductivo Básico: $`R_0 = \frac{\beta}{\gamma} > 1`$.En la gráfica,
para $`\beta=0.1`$ (donde $`\beta = \gamma`$), la curva es plana o
decrece (no hay epidemia). A medida que $`\beta`$ aumenta, el pico es
más alto y ocurre más pronto.

\##Pregunta 4 Variación de $`\gamma`$ (Recuperación) Estudiamos la
dinámica manteniendo $`\beta = 1`$ y variando $`\gamma`$.Valores de
$`\gamma`$: 0.025, 0.2, 0.5, 1.

``` r
beta_fijo <- 1
gammas <- c(0.025, 0.2, 0.5, 1)
resultados_p4 <- list()

for(g in gammas){
  pars <- c(beta = beta_fijo, gamma = g)
  out <- as.data.frame(ode(y=initial_state_values, times=times, func=sir_model, parms=pars))
  out$gamma <- as.factor(g)
  resultados_p4[[length(resultados_p4)+1]] <- out
}

data_p4 <- do.call(rbind, resultados_p4)

ggplot(data_p4, aes(x=time, y=I, color=gamma)) +
  geom_line(size=1) +
  labs(title = "Dinámica variando Gamma (Beta = 1)", 
       y = "Infectados", color = "Gamma") +
  theme_minimal()
```

![](Reto_Parte_1_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
*Análisis*: Al aumentar $`\gamma`$ (recuperación más rápida), el pico de
infectados disminuye.Si $`\gamma = 1`$, entonces $`\beta = \gamma`$ y
$`R_0 = 1`$. En este caso (línea correspondiente a 1), no se observa un
brote epidémico significativo.Nuevamente se confirma que para que exista
una epidemia, la tasa de recuperación debe ser suficientemente baja
comparada con la de infección: $`\gamma < \beta`$.
