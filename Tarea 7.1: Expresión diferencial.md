## TAREA 7.1: Expresión diferencial
**Curso:** Introducción a la bioinformática e investigación reproducible para análisis genómicos

**Unidad 6:** Introducción a la genómica y secuenciación de siguiente generación

**Catedráticos:** Dr. Ricardo Verdugo-Salgado y Dra. Karen Oróstica

**Estudiante:** José Daniel Lara-Tufiño

**Reporte de tarea 2.1:** Expresión diferencial

**Fecha de entrega:** 29/04/20

### INTRODUCCIÓN

La detección de genes expresados diferencialmente bajo distintas condiciones experimentales es un proceso crucial en cualquier anális de secuenciación de nueva generación (NGS). La efectividad de las mediciones de expresión genética está altamente correlacionada con el diseño experimental. En este sentido, la aleatorización es reconocida universalmente como un elemento indispensable en cualquier diseño experimental, pues permite evitar los sesgos de selección (a distintos tratamientos por ejemplo). Sin embargo, la aleatorización a veces se ignora en algunas etapas de la obtención, procesesamiento y análisis de los datos de NGS; éste es el caso en los diseños experimentales con *microarrays*, en los cuales existen varias alternativas para colocar las muestras en los microarreglos, colocandose algunas veces agrupadas, por ejemplol: controles vs unidades con tratamiento, por sexo, especie entre otras.

El estudio del cual se obtuvieron los datos, evaluó dos experimentos, el primero colocó el genotipo con el BeadChip y el tratamiento con la posición de la matriz, mientras que en el segundo éstas variables se asignaron al azar con el BeadChip y la posición de la matriz. Finalmente, en el presente reporte se evaluó el efecto de la variación genética en el cromosoma Y del ratón sobre el tamaño de los cardiomiocitos y la posible dependencia de tales efectos en niveles de testosterona.

### OBJETIVOS

**1-** Determinar si existe expresión diferencial entre genotipos.

**2-** Determinar si existe expresión diferencial entre tratamientos.

**3-** Evaluar las diferencias en la respuesta al tratamiento entre los dos genotipos.

### MÉTODOS: Diseño experimental

Se evaluaron ocho ratones machos adultos de dos cepas (C57BL/6J y C57BL/6J-chrY<A/J/NaJ>9), denominadas **B** y **BY** de aquí en adelante. De cada cepa (genotipo), cuatro animales fueron castrados y cuatro fueron intervenidos con el mismo procedimiento quirúrjuico, excepto que no se realizó la castración (animales intactos usados como control). El ARN se hibridizó a BeadChips Illumina MouseRef-8 v2.0 que contienen ocho microarreglos con 25,697 sondas cada uno. Unicamente se seleccionaron 5000 sondas de manera aleatoria para este tutorial.

## Indicaciones: 

1- Obtener la matriz completa de datos desde GEO, ver link en la Descripción de los datos e importar la matriz en R, seleccionar aleatoriamente 5000 filas y exportar el subset de datos en un archivo plano separado por tabulaciones.

2- Ejecute este tutorial, pero con algunos cambios:

-Usando su matriz de datos en vez de la usada en la demostración.
-En vez de considerar un transcrito presente si la sonda lo detectó en el 50% de las muestras de cualquier grupo experimental, hágalo cuando se detectó en al menos el 25% de las muestras de todos grupos experimentales.
-En vez de usar 200 permutaciones, use 500.
-En vez de usar un FDR 0.2, use uno de 0.19.
-En vez de seleccionar un gen si cualquier sonda asociada al gen está seleccionada, hágalo solo cuando todas las sondas lo están.

**Nota**: No se muestra una serie de comandos exactamente secuencial, sólo se reportan los comandos y resultados que se consideraron más importantes.

### Descargar 500 muestras de manera aleatoria del archivo GSE15354_raw.txt obtenido de http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15354. También, leer los valores crudos (no normalizados):

```R
Data <- read.delim("C:/Raiz/BioinfinvRepro-master/Unidad7/DE_tutorial/raw_data.txt", header=TRUE)
my_data <- Data[sample(nrow(Data), 5000), ]
signal    <- grep("CDR", colnames(my_data)) # vector de columnas con datos 
detection <- grep("Detection.Pval", colnames(my_data)) # vector de columnas con valores p
```

### Importar las anotaciones de las sondas y extraer las 5000 previamente seleccionadas:

```R
annot     <- read.delim("C:/Raiz/BioinfinvRepro-master/Unidad7/DE_tutorial/MouseRef-8_annot_full.txt")
data_names <- rownames(my_data)
annot = annot[row.names(annot)%in%data_names,]
```

### Visualizar la calidad de las sondas al ser alineadas con el genoma de referencia: 

```R
table(annot$ProbeQuality)

Bad        Good      Good***    Good****    No match     Perfect     Perfect***     Perfect****
273         101       4           14            7          4429         51             121
```

### Tabla con el diseño de hibridaciones:

```R
design <- read.csv("C:/Raiz/BioinfinvRepro-master/Unidad7/DE_tutorial/YChrom_design.csv")
print(design)

Array Sample_Name Sentrix_ID Sentrix_Position Genotype Treatment Group
1      1  CDR017-DIL 4340571022                A        B         I   B.I
2      2      CDR025 4340571022                B       BY         I  BY.I
3      3      CDR021 4340571022                C        B         C   B.C
4      4  CDR029-DIL 4340571022                D       BY         C  BY.C
5      5      CDR022 4340571022                E        B         C   B.C
6      6      CDR018 4340571022                F        B         I   B.I
7      7  CDR026-DIL 4340571022                G       BY         I  BY.I
8      8  CDR030-DIL 4340571022                H       BY         C  BY.C
9      9  CDR031-DIL 4340571033                A       BY         C  BY.C
10    10      CDR023 4340571033                B        B         C   B.C
11    11  CDR027-DIL 4340571033                C       BY         I  BY.I
12    12      CDR019 4340571033                D        B         I   B.I
13    13      CDR020 4340571033                E        B         I   B.I
14    14  CDR028-DIL 4340571033                F       BY         I  BY.I
15    15  CDR032-DIL 4340571033                G       BY         C  BY.C
16    16      CDR024 4340571033                H        B         C   B.C   
```
## Control de calidad
### Crea gráficos de cajas coloreados por la calidad de la sonda:

```R
png(file.path(outdir,"boxplot_raw_probe_qc.png"), width=6.5, height=4, unit="in", res=150)

par(xpd=NA, mar= c(6.1, 5.1, 4.1, 2.1), cex=.7, las=3)
boxplot(unlist(log2(my_data[,signal]))~qcfact+afact, horiz=T, main="Raw log2 values Boxplot",
        col=rep(1:2, length(signal)), axes=F, varwidth=TRUE)
axis(1, at=seq(1, length(signal)*2, by=2)+.5, labels=alabel)
axis(2)
legend("top", legend=levels(qcfact), fill=1:2, ncol=2, xjust=.5, bty="n", inset=-.1)

dev.off()
```
![bloxplot_raw_probe_qc](https://github.com/jdaniellt/Tarea-7.1-Expresi-n-diferencial/blob/master/boxplot_raw_probe_qc.png)

### Crea cuadros de caja de colores por tratamiento

```R
png(file.path(outdir,"boxplot_raw_treatment.png"), width=4, height=4, unit="in", res=150)

par(xpd=NA, mar= c(6.1, 4.1, 4.1, 2.1), cex=.7)
boxplot(as.data.frame(log2(my_data[,signal])), horiz=T, main="Raw log2 values Boxplot", las=1, col=design$Treatment, names=design$Sentrix_Position, cex.axis=.9)
legend(8, 2.5, legend=levels(design$Treatment), fill=1:2, ncol=2, xjust=.5)

dev.off()
```
![boxplot_raw_treatment](https://github.com/jdaniellt/Tarea-7.1-Expresi-n-diferencial/blob/master/boxplot_raw_treatment.png)




