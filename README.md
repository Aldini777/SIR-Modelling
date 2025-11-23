#  Simulación de la Evolución de Contagios (Modelo SIR)

> **Materia:** Sistemas de Ecuaciones Diferenciales  
> **Autor:** Aldo Resendiz  
> **Herramientas:** R, RStudio / VS Code

---

##  Introducción

Este repositorio contiene el desarrollo del Reto final para la materia de **Sistemas de Ecuaciones Diferenciales**. El objetivo principal es la modelación matemática y computacional de fenómenos epidemiológicos para entender cómo se propagan las enfermedades en poblaciones controladas.

##  Evolución de Contagios

En el reto de este bloque, consideraremos el fenómeno de evolución de contagios de una enfermedad en un sistema poblacional cerrado, donde las posiciones y movimientos de los individuos tienen una naturaleza aleatoria. 

Este tipo de modelos, llamados **SIR** (Susceptibles - Infectados - Recuperados), pueden servir para:
* Pronosticar las curvas de contagios de una enfermedad.
* Tomar medidas para corregir la evolución de manera que se tengan condiciones adecuadas.

Es precisamente en estos tiempos de pandemia en los que la utilidad de este tipo de modelos se vuelve relevante a nivel global.

### Enfoque del Proyecto
Los contenidos de los módulos referentes a las consideraciones numéricas y teóricas deberán ser utilizadas y complementadas con consideraciones de **tipo estocástico** para determinar la dinámica de contagios en un ambiente controlado.

La situación a considerar incluye un grupo de población de una especie que es afectada por una enfermedad infecciosa y cuyo comportamiento dinámico de movimiento no puede ser conocido con antelación, puesto que cada individuo puede moverse siguiendo algunas reglas básicas de movimientos posibles.

Con esta dinámica desconocida, será necesario plantear y resolver un **sistema de ecuaciones diferenciales** para describir a los:
1. Individuos **Enfermos**.
2. Individuos **Susceptibles** de contraer la enfermedad.
3. Individuos **Recuperados** o removidos de la dinámica contagiosa.

---

##  Etapas del Proyecto

Con el fin de generar una solución de complejidad creciente, este proyecto se divide en cuatro etapas fundamentales:

1.  **Modelación y simulación básica:** Difusión de enfermedades con el modelo SIR estándar.
2.  **Refinamiento del modelo:** Modificaciones del modelo SIR para considerar mejoras (dinámica vital, vacunación, etc.).
3.  **Dinámica Estocástica:** Consideración de la aleatoriedad en la simulación de contagios.
4.  **Simulación Espacio-Temporal:** Simulación de la dinámica de los contagios considerando el movimiento de los individuos del sistema.

---

##  Visualización de Resultados

Para ver los gráficos y el análisis del código renderizado directamente en GitHub, por favor abre los archivos con extensión `.md` (por ejemplo, `Reto.md`) en este repositorio.
