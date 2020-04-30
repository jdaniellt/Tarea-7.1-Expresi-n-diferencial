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

### Cargar las siguientes librerias:

```R 
library(preprocessCore)
library(maanova)
library(limma)
library(topGO)
library(org.Mm.eg.db)
```

### Defina algunas constantes. En vez de usar un FDR 0.2, use uno de 0.19:

```R 
outdir     <- "output"
fdr_th     <- 0.19
```

### Lea un archivo que define algunas funciones que son necesarias para el análisis:

```R
> source("Rfxs.R")
```

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
### Figura 1:
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
### Figura 2
```R
png(file.path(outdir,"boxplot_raw_treatment.png"), width=4, height=4, unit="in", res=150)

par(xpd=NA, mar= c(6.1, 4.1, 4.1, 2.1), cex=.7)
boxplot(as.data.frame(log2(my_data[,signal])), horiz=T, main="Raw log2 values Boxplot", las=1, col=design$Treatment, names=design$Sentrix_Position, cex.axis=.9)
legend(8, 2.5, legend=levels(design$Treatment), fill=1:2, ncol=2, xjust=.5)

dev.off()
```
![boxplot_raw_treatment](https://github.com/jdaniellt/Tarea-7.1-Expresi-n-diferencial/blob/master/boxplot_raw_treatment.png)

### Cree diagramas de dispersión de datos sin procesar en escala log2.
### Figura 3
```R
png(file.path(outdir,"Pairs_scatter_log2.png"), width=8, height=8, unit="in", res=150)
par(cex=.2, mar=c(2.1,2.1,2.1,1.1))
pairs(log2(Data.Raw[,signal]), main="Log2 Raw Intensity Values", pch=".",  gap=.5, cex.labels=.5)
dev.off()
```
![Pairs_scatter_log2](https://github.com/jdaniellt/Tarea-7.1-Expresi-n-diferencial/blob/master/Pairs_scatter_log2.png)

### Filtrado de la sonda por QC: 
**Las sondas de mala calidad tienden a tener una señal más baja que las sondas buenas. Por lo tanto se recomienda para eliminarlas**

```R
my_data <- Data.Raw[probe_qc %in% "Good probes",]
annot    <- annot[probe_qc %in% "Good probes",]
```
### Crear un microarreglo de datos brutos

```R 
rawdata           <- as.matrix(my_data[,signal])
rownames(rawdata) <- Data.Raw$PROBE_ID
colnames(rawdata) <- design$Sample_Name
```

### Normalización de datos:

```R
library(preprocessCore)
normdata           <- normalize.quantiles(rawdata) 
colnames(normdata) <- colnames(rawdata)
rownames(normdata) <- rownames(rawdata)
```

### Filtrado de sondas. Este paso tiene como objetivo eliminar las sondas que no detectaron trasncritos en ninguno de los grupos experimentales. Tenga en cuenta que este paso puede ser opcional pero recomendado. 

**En vez de considerar un transcrito presente si la sonda lo detectó en el 50% de las muestras de cualquier grupo experimental, hágalo cuando se detectó en al menos el 25% de las muestras de todos grupos experimentales:** 

```R 
probe_present      <- my_data[,detection] < 0.04
detected_per_group <- t(apply(probe_present, 1, tapply, design$Group, sum))
present  <- apply(detected_per_group >= 1, 1, all)
normdata <- normdata[present,]
annot    <- annot[present, ]
```
## Prueba de expresión diferencial

```R
 madata <- read.madata(normdata, design, log.trans=T)
```
### Ajustar el modelo. 
```R
 fit.fix <- fitmaanova(madata, formula=~Group)
```
### Estime algunas estadísticos básicos para cada grupo experimental para ser incluidas en la tabla final de los resultados. 
```R
 Means           <- t(apply(madata$data, 1, tapply, design$Group, mean)) 
 colnames(Means) <- paste("Mean", colnames(Means), sep=":")
 SEs             <- t(apply(madata$data, 1, tapply, design$Group, function(x) sqrt(var(x)/length(x))))
 colnames(SEs)   <- paste("SE", colnames(SEs), sep=":")
```
### Pruebe cada contraste utilizando 500 permutaciones de las muestras. En una situación real se recomiendan al menos 1.000 permutaciones. Las pruebas de F se realizarán utilizando una estimación de varianza residual por sonda (F1) y una estimación basada en contracción de varianza residual que utiliza información de múltiples sondas (Fs)

```R 
test.cmat <- matest(madata, fit.fix, term="Group", Contrast=cmat, n.perm=500, 
                    test.type = "ttest", shuffle.method="sample", verbose=TRUE)
                    Doing F-test on observed data ...
Doing permutation. This may take a long time ... 
Finish permutation #  100 
Finish permutation #  200 
```
### Grafique los valores de p comparando diferentes formas de calcularlos ( consulte ?matest):
### Figura 4
```R
png(file.path(outdir,"P-values Hist.png"), width=6, height=6, unit="in", res=150)
par(mfrow=c(2,2), oma=c(2,0,2,0), cex=.8, xpd=NA)
palette(rainbow(3))
plot(density(test.cmat$F1$Ptab[,1]), col=1, main="F1:Ptab", lwd=2)
lines(density(test.cmat$F1$Ptab[,2]), col=2, lwd=2)
lines(density(test.cmat$F1$Ptab[,3]), col=3, lwd=2)

plot(density(test.cmat$F1$Pvalperm[,1]), col=1, main="F1:Pvalperm", lwd=2)
lines(density(test.cmat$F1$Pvalperm[,2]), col=2, lwd=2)
lines(density(test.cmat$F1$Pvalperm[,3]), col=3, lwd=2)

plot(density(test.cmat$Fs$Ptab[,1]), col=1, main="Fs:Ptab", lwd=2)
lines(density(test.cmat$Fs$Ptab[,2]), col=2, lwd=2)
lines(density(test.cmat$Fs$Ptab[,3]), col=3, lwd=2)

plot(density(test.cmat$Fs$Pvalperm[,1]), col=1, main="Fs:Pvalperm", lwd=2)
lines(density(test.cmat$Fs$Pvalperm[,2]), col=2, lwd=2)
lines(density(test.cmat$Fs$Pvalperm[,3]), col=3, lwd=2)

legend(-.5, -1.6, legend=c("Geno", "Trt", "Int"), col=1:3,lwd=2,xjust=.5,ncol=3,xpd=NA)
dev.off()
```

![P-values%20Hist](https://github.com/jdaniellt/Tarea-7.1-Expresi-n-diferencial/blob/master/P-values%20Hist.png)

###  Resuma en una tabla los resultados para todas los transcriptos presentes.  Exportaremos solo los resultados de las pruebas `Fs` y los valores `Pvalperm` de permutaciones.

```R
 results <- data.frame(annot, Means, SEs, F_val=test.cmat$Fs$Fobs,
                       P_val=test.cmat$Fs$Pvalperm, FDR=test.cmat$Fs$adjPvalperm, FC=FC)
```

### Exportar todos los resultados.  Puede abrir DE_results.csv en Excel o Calc OpenOffice.

```R
 write.table(results, file=file.path(outdir,"DE_results.csv"), sep=",", row.names=F)
```

## Contar genes expresados diferencialmente

En este experimento, una pregunta de interés fue cuántos genes responden de manera diferente al tratamiento de la castración en los dos genotipos. En otras palabras, ¿existe un efecto de interacción entre el genotipo y el tratamiento sobre la expresión génica en los cardiomiocitos? En segundo lugar, es interesante evaluar la naturaleza de la interacción. ¿Ambos genotipos responden al tratamiento pero en direcciones opuestas? ¿O el tratamiento tiene un efecto en uno de los genotipos y no en el otro? Para responder a la primera pregunta, se podría contar el número de sondas que muestran un efecto significativo para el contraste Int, es decir, el número de sondas que tienen un FDR por debajo de un umbral para la prueba Int.Pvalperm. Si aún no sabe cómo calcular esto en R, abra el archivo DE_results.csv un editor de hoja de cálculo como Calc OpenOffice y use filtros para calcular este número. Luego trate de hacer esto en R y compare los resultados.

Para responder la segunda pregunta, podemos usar los diagramas de Venn. Queremos contar el número de genes que se seleccionan en cada genotipo, pero solo entre aquellos que muestran un efecto de interacción significativo.

### Crear un identificador del gen basado en EntrezGene y utilizar el ID de la sonda cuando no esté asociada a un gen:

```R
 results$GeneID <- results$EntrezID
 results$GeneID[is.na(results$GeneID)] <- results$ProbeID[is.na(results$GeneID)]
 ```
 
### Cuente las sondas seleccionadas por expresión diferencial por genotipo, tratamiento y/o interacción:

```R
 Probes.DE <- results[, c("FDR.Geno", "FDR.Trt", "FDR.Int")]  <= fdr_th
 Genes.DE  <- apply(Probes.DE, 2, tapply, results$GeneID, any)
```

### Usando solo las sondas seleccionadas por efectos de interacción, cuente las sondas significativas para el efecto de genotipo en ratones intactos (I) y/o castrados (C).

```R             
 Probes.Int_Geno <- results[results$FDR.Int <= fdr_th, 
                            c("FDR.Geno_I", "FDR.Geno_C")] <= fdr_th
 Genes.Int_Geno  <- apply(Probes.Int_Geno, 2, tapply, 
                          results$GeneID[results$FDR.Int <= fdr_th], any)
   ```
### Usando solamente sondas seleccionadas por efectos de interacción, cuente las sondas significativas para el efecto de tratamiento en ratones del genotipo B y/o del genotipo BY.

```R
 Probes.Int_Trt  <- results[results$FDR.Int <= fdr_th,
                            c("FDR.Trt_B", "FDR.Trt_BY")]  <= fdr_th
 Genes.Int_Trt   <- apply(Probes.Int_Trt, 2, tapply,
                          results$GeneID[results$FDR.Int <= fdr_th], any)
```
                          
### Contar genes para cada combinación de efectos marginales y de interacción.

```R 
Counts.DE <- vennCounts(Genes.DE)
print(Counts.DE)
  FDR.Geno FDR.Trt FDR.Int Counts
1        0       0       0    863
2        0       0       1     36
3        0       1       0    437
4        0       1       1     33
5        1       0       0    467
6        1       0       1     30
7        1       1       0    302
8        1       1       1     56
attr(,"class")
[1] "VennCounts"
```

### Contar los genes DE entre niveles de un factor condicional en el otro factor.

```R
Counts.Int_Geno <- vennCounts(Genes.Int_Geno)
print(Counts.Int_Geno)
  FDR.Geno_I FDR.Geno_C Counts
1          0          0      5
2          0          1     97
3          1          0     39
4          1          1     14
attr(,"class")
[1] "VennCounts"
```
```R 
Counts.Int_Trt  <- vennCounts(Genes.Int_Trt) 
> 
> print(Counts.Int_Trt)
  FDR.Trt_B FDR.Trt_BY Counts
1         0          0      2
2         0          1    109
3         1          0     27
4         1          1     17
attr(,"class")
[1] "VennCounts"
```

### Graficar los genes DE por efectos marginales o de interacción. 
### Figura 5
```R
png(file.path(outdir, "vennDiagram_DiffExprs.png"), width=3.5, height=3, unit="in", res=150)
 par(cex=.7)
 vennDiagram(Counts.DE, names=c("Geno", "Trt", "Int"), 
             main="\n\n\nDifferentially Expressed Genes")
 dev.off()
 ```

![vennDiagram_DiffExprs](https://github.com/jdaniellt/Tarea-7.1-Expresi-n-diferencial/blob/master/vennDiagram_DiffExprs.png)

### Generación de Diagramas de Venn de los genes que responden al genotipo de manera dependiente del tratamiento y viceversa
### Figura 6
```R
png(file.path(outdir, "vennDiagram_Int.png"), width=6.5, height=3, unit="in", res=150)

 par(mfrow=c(1,2), cex=.7)
 vennDiagram(Counts.Int_Geno, names=c("I", "C"), 
             main="\n\n\nGenes Responding to Genotype\nin a Treatment Dependent Manner")
 vennDiagram(Counts.Int_Trt, names=c("B", "BY"),
             main="\n\n\nGenes Responding to Treatment\nin a Genotype Dependent Manner")

 dev.off()
```

![vennDiagram_Int](https://github.com/jdaniellt/Tarea-7.1-Expresi-n-diferencial/blob/master/vennDiagram_Int.png)


## Pruebas funcionales

### Las sondas no anotadas no se pueden usar para el análisis de enriquecimiento, así que las eliminaremos

```R 
results   <- results[!is.na(results$EntrezID),]
```

### Sondas seleccionadas por interacción

```R 
probes.int <- results[, "FDR.Int"]  <= fdr_th
```

### Resumir por gen. Se selecciona un gen si se selecciona cualquier transcripto.

```R 
genes.int <- tapply(probes.int, results$EntrezID, any)
```

### Convertir de lógico a vector de 0 (no seleccionado) o 1 (seleccionado)

```R
genes.int <- ifelse(genes.int, 1, 0)
```

### Convertir a un factor

```R
genes.int <- as.factor(genes.int)
```
### Cree un objeto de datos de clase topGO, que contiene todo lo necesario para el enriquecimiento de las pruebas.

```R 
 GOdata <- new("topGOdata", ontology="BP", allGenes=genes.int,
               description="Genes DE by Trt by GenoInteraction", nodeSize=5,
               annotationFun=annFUN.org, mapping="org.Mm.eg.db", ID="entrez")  
```

### Pruebe si existe enriquecimiento de términos GO utilizando una prueba exacta de Fisher. Usaremos dos algoritmos, la prueba clásica de término por término y el algoritmo " elim " , que tiene en cuenta la jerarquía de los términos GO para evitar la redundancia.

```R 
resultFisher.classic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
 resultFisher.elim    <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

 GO_BP_Table <- GenTable(GOdata, topNodes = 20, Fisher.classic=resultFisher.classic, 
                         Fisher.elim=resultFisher.elim, 
                         orderBy = "Fisher.elim", ranksOf = "Fisher.classic")
                         
  GO.ID                                        Term Annotated Significant
1  GO:0032456                         endocytic recycling         5           4
2  GO:0036010            protein localization to endosome         5           4
3  GO:0043393               regulation of protein binding        35           8
4  GO:0051963              regulation of synapse assembly        16           5
5  GO:0061951 establishment of protein localization to...        10           4
6  GO:0098876 vesicle-mediated transport to the plasma...        11           4
7  GO:1903649         regulation of cytoplasmic transport         6           3
8  GO:1903670        regulation of sprouting angiogenesis         6           3
9  GO:1905477 positive regulation of protein localizat...        13           4
10 GO:0035306    positive regulation of dephosphorylation         7           3
11 GO:0043537 negative regulation of blood vessel endo...         7           3
12 GO:0016197                         endosomal transport        26           9
13 GO:0014706          striated muscle tissue development        57           9
14 GO:0010811 positive regulation of cell-substrate ad...        21           5
15 GO:0110020 regulation of actomyosin structure organ...        14           4
16 GO:0051495 positive regulation of cytoskeleton orga...        30           6
17 GO:0048468                            cell development       289          28
18 GO:0032091      negative regulation of protein binding        15           4
19 GO:0055007         cardiac muscle cell differentiation        23           5
20 GO:0051146        striated muscle cell differentiation        41           7
   Expected Rank in Fisher.classic Fisher.classic Fisher.elim
1      0.32                      2        7.4e-05     7.4e-05
2      0.32                      3        7.4e-05     7.4e-05
3      2.22                      5         0.0012      0.0012
4      1.02                      6         0.0024      0.0024
5      0.64                      7         0.0024      0.0024
6      0.70                      8         0.0036      0.0036
7      0.38                      9         0.0044      0.0044
8      0.38                     10         0.0044      0.0044
9      0.83                     13         0.0071      0.0071
10     0.44                     14         0.0073      0.0073
11     0.44                     15         0.0073      0.0073
12     1.65                      1        1.6e-05      0.0076
13     3.62                     18         0.0085      0.0085
14     1.33                     19         0.0086      0.0086
15     0.89                     20         0.0095      0.0095
16     1.91                     21         0.0098      0.0098
17    18.36                     25         0.0116      0.0116
18     0.95                     26         0.0123      0.0123
19     1.46                     27         0.0128      0.0128
20     2.61                     28         0.0130      0.0130
```
### Exporte la tabla completa de los principales términos del formato CSV de GO (puede abrirse en Excel o Calc OpenOffice ) .

```R
 write.table(GO_BP_Table, file.path(outdir, "GO_BP_Table.csv"), sep=",", 
  row.names=F)
```

## RESULTADOS Y DISCUSIÓN
Existe un mayor número de genes DE entre genotipos (527) que entre tratamientos (428), aunque una considerable cantidad de genes DE son compartidos (309) en ambos casos (Fig. 5). Con referencia a las interacciones, se observa un resultado similar, hay más genes DE entre los genotipos (37) que en entre tratamientos (30), y se comparten 53 genes DE en ambos casos (Fig. 5). Como producto de las interacciones, se identificaron menos genes DE (23), con relación a los factores experimentales (genotipos y tratamientos).Con Referencia a los genes que responden al genotipo de manera dependiente del tratamiento y viceversa, se obtuvieron resultados más contrastantes (Fig. 6). Se identificó que los ratones intactos presentan un menor número de genes ED (29) que los ratones castrados (88), y con referencia al genotipo, la cepa BY muestra un número mayor de genes DE (22) que la cepa B (92).

En cuanto a las implicaciones biológicas, la evaluación de la interacción permitió identificar 249 “genes de interacción”  que se enriquecen dentro la vía de señalización p53, que comprende los genes Ccnd1, Pten y Cdkn1a que son correguladores de los receptores de andrógenos, así como la vía de señalización del ritmo circadiano que comprende los genes Arntl / Bmal1, encargados de regular el crecimiento celular a través del control de Cdkn1a.

Se encontró, que la testosterona post-puberal aumentó el tamaño de los cardiomiocitos de la cepa C56BL / 6J, pero no fue así con la cepa C57BL / 6J-chrYA; sin embargo, no afectó la expresión génica en los corazones de ambas. No obstante, varios genes cardíacos respondieron a la testosterona post-puberal de una manera estrictamente selectiva dependiendo de la cepa. Lo anterior, puede explicarse por la ausencia de efectos hipertróficos de post-pubertad de la testosterona en cada una de las cepas. Finalmente, los hallazgos brindan evidencia de que los genes de los cromosomas sexuales juegan papel importante en las enfermedades cardiovasculares.











