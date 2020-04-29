## TAREA 7.1: Expresión diferencial
**Curso:** Introducción a la bioinformática e investigación reproducible para análisis genómicos

**Unidad 6:** Introducción a la genómica y secuenciación de siguiente generación

**Catedráticos:** Dr. Ricardo Verdugo-Salgado y Dra. Karen Oróstica

**Estudiante:** José Daniel Lara-Tufiño

**Reporte de tarea 2.1:** Expresión diferencial
m
**Fecha de entrega:** 29/04/20

### INTRODUCCIÓN

La detección de genes diferencialmente expresados bajo diferentes condiciones experimentales, es un proceso crucial en cualquier anális de secuenciación de nueva generación.
observar si existen distintos comportamientos en la producción de proteínas en las células cuando
éstas se encuentran sometidas a estímulos. Estos estímulos pueden ser de dos tipos: internos,
como los que derivan en la diferenciación celular, o externos debidos al entorno. Muchas de
las herramientas más utilizadas para llevar a cabo este análisis realizan comparaciones dos a
dos entre las condiciones experimentales (típicamente, condición de control frente a condición
experimental) [4–6] y no permiten identificar de manera sencilla aquellos genes que se expresan
diferencialmente en una sola condición experimental [7]. Sin embargo, cada día es más común
tratar de identificar estos genes específicos que nos pueden indicar por ejemplo alteraciones
celulares debidas a algún cáncer [8] o incluso la presencia de diferentes toxinas en el ambiente.
Actualmente hay pocas técnicas destinadas a la detección de genes diferencialmente expresados bajo varias condiciones experimentales [8–12], y además ninguna se centra en la especificidad
de ese gen con respecto al resto de condiciones, por lo que es un problema abierto dentro del
campo de la bioinformática que precisa de soluciones eficaces. Además el abaratamiento del
la tecnología de secuenciación de la última década está dando lugar a estudios experimentales
cada vez más complejos y que precisan soluciones más específicas. Un ejemplo de aplicación es
la detección e identificación de diferentes sustancias tóxicas en el agua mediante la bacteria E.
Coli [13].

### OBJETIVOS

**1-**Determinar si existe expresión diferencial entre genotipos.

**2-**Determinar si existe expresión diferencial entre tratamientos.

**3-**Evaluar las diferencias en la respuesta al tratamiento entre los dos genotipos.

### MÉTODOS: Diseño experimental

Se evsluaron ocho ratones machos adultos de dos cepas (C57BL/6J y C57BL/6J-chrY<A/J/NaJ>9), denominadas **B** y **BY** de aquí en adelante. De cada cepa (genotipo), cuatro animales fueron castrados y cuatro fueron intervenidos con el mismo procedimiento quirúrjuico, excepto que no se realizó la castración (animales intactos usados como control). El ARN se hibridizó a BeadChips Illumina MouseRef-8 v2.0 que contienen ocho microarreglos con 25,697 sondas cada uno. Unicamente se seleccionaron 5000 sondas de manera aleatoria para este tutorial.






Figura 1: 
![alt text](https://github.com/jdaniellt/Figura-1.-Diagrama-de-caja-de-datos-sin-procesar-en-escala-log-por-microarreglo-y-calidad-de-sonda./blob/master/Pairs_scatter_log2.png)
Raquelita no es enojoncita

Es suavecita
