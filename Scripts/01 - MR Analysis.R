
# 00 - Install and load packages ------------------------------------------

# Packages 
#install.packages("MendelianRandomization", dependencies = TRUE)
#install.packages("metafor", dependencies = TRUE)
library(MendelianRandomization)
library(metafor)

##Data that comes with the MendelianRandomization package
#usethis::use_data(calcium)
#usethis::use_data(ldlc)

# install MRbase items; will use this for pruning SNPs in LD
# install.packages("devtools")
# library(devtools)
# install_github("MRCIEU/TwoSampleMR")
# library(TwoSampleMR)

#install MRPRESSO; identifies IVs with heterogeneous SNPs
#if (!require("devtools")) { install.packages("devtools") } else {}
#devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)

#install.packages("tidyverse", dependencies = TRUE)
library(tidyverse)


# 01 - GWAS of uric acid in Hispanic Community Health Study/ Study of Latinos (HCHS/SOL) ----

# Exposure GWAS “SOL_urate_summary.csv.gz”  
# GWAS of uric acid (Feofanova et al AJHG 2020) 

sol_urate <- read_csv(gzfile("Data/SOL_urate_summary.csv.gz"), col_name=TRUE)

## 1.1 
# Seleccionamos las variables con un valor de P < 5x10^{-6}
sol_urate_lim <- sol_urate[sol_urate$PVAL < 0.000005, ] 
nrow(sol_urate_lim) # 1169 variantes 

## 1.2 
# Creamos una variable de identificación ID con el cromosoma, posición, y alelos
# para unir con la base de datos de evento

sol_urate_lim <- 
  sol_urate_lim %>%
  mutate(
    cpaid = str_c(chromosome,                 # chromosome 
                  position,                   # pos
                  pmin(alleleA, alleleB),     # alelo menor
                  pmax(alleleA, alleleB),     # alelo mayor
                  sep=":")                    # Separador
    )         

## 1.3 
# Renombramos las columnas en el archivo de exposición para homogeneizar las
# bases de datos

exposure <-
  data.frame(cpaid = sol_urate_lim$cpaid,
             beta.exposure = sol_urate_lim$Beta,
             se.exposure = sol_urate_lim$SE,
             effect_allele.exposure = sol_urate_lim$alleleA,
             other_allele.exposure =sol_urate_lim$alleleB,
             eaf.exposure = sol_urate_lim$AF,
             pval.exposure = sol_urate_lim$PVAL)

# Se mantienen la mayor parte de las bases de datos con expeción de la pos y el cromosoma

head(exposure)

# 02 -  Cardiogram consotrium, GWAS of Coronary artery disease ------------

# Outcome GWAS “cad.add.160614.website.txt.gz”
# Nikpay et al, CAD Nature Genet, 2015

cad <- read.table(gzfile("Data/cad.add.160614.website.txt.gz"), header = TRUE)

nrow(cad)

## 2.1
# Se genera un ID para las variantes de urato con cromosoma, posición, y
# alelos para unir con la base de datos de exposición

cad <- 
  cad %>% 
  mutate(
    cpaid = str_c(chr,                                      # chromosome 
                  bp_hg19,                                  # pos
                  pmin(effect_allele, noneffect_allele),    # alelo menor
                  pmax(effect_allele, noneffect_allele),    # alelo mayor
                  sep=":")                                  # Separador
  )

## 2.3 
# Seleccionamos el instrumento genético SNP's de la base de datos de evento 

sl <- cad[cad$cpaid %in% exposure$cpaid, ] # get variants in CAD GWAS

nrow(sl)

## 2.4 
# Renombramos las columnas en la base de datos de evento para que sean los
# mismos en todos los archivos

outcome <- 
  data.frame(cpaid = sl$cpaid,
             beta.outcome = sl$beta,
             se.outcome = sl$se_dgc,
             effect_allele.outcome = sl$effect_allele,
             other_allele.outcome = sl$noneffect_allele,
             eaf.outcome = sl$effect_allele_freq,
             pval.outcome = sl$p_dgc,
             SNP = sl$markername)

head(outcome)


# 03 - Unión de las bases de datos ----------------------------------------

# Unimos las bases de datos de exposición y yde evento  
exp_out <- 
  merge(outcome, exposure, 
        by = "cpaid",
        all.x = TRUE) 

nrow(exp_out)
head(exp_out)

# 04 - SNPs independientes ------------------------------------------------

# Usamos MRbase para encontrar SNPs independientes con r2 < 0.05 (TwoSampleMR) 

all <- TwoSampleMR::clump_data(exp_out, 
                  clump_r2 = 0.05, 
                  pop = "AMR"             # Population reference panel - America 
                  ) 

nrow(all)

# Quitamos el paquete TwoSampleMR de nuestro workflow ya que impide el trabajo
# con al paquete MendelianRandomization

# detach(package:TwoSampleMR)


# Checkpoint --------------------------------------------------------------

# save.image(file = "Output/202309 - MR Analysis_RData.RData")

# Antes de continuar establecemos un punto de control desde el cual se puede 
# cargar todas las bases de datos anteriores en el área de trabajo

# load("Output/202309 - MR Analysis_RData.RData")


# 05 - Estadística F ------------------------------------------------------

all$f1 <-(all$beta.exposure * all$beta.exposure) / (all$se.exposure * all$se.exposure)
all$f1

mean(all$f1)

# 06 - Alineación de los SNPs con mismo alelo de efecto -------------------

# Alineamos los SNPs con el mismo alelo de efecto para la exposición y el evento
# al cambiar el el símbolo de la beta de la exposición si los efectos de los
# alelos no son similares.

# Convertimos los valores de los alelos de exposición y de evento a factores
all$effect_allele.outcome <- as.factor(all$effect_allele.outcome)
all$effect_allele.exposure <- as.factor(all$effect_allele.exposure)

# Se extren todos los alelos de exposición y de evento evitando valores repetidos 
lev2 <- unique(c(levels(all$effect_allele.outcome), levels(all$effect_allele.exposure)))

# Asignamos los niveles a cada variable  
all$effect_allele.outcome <- factor(all$effect_allele.outcome, levels=lev2)
all$effect_allele.exposure <- factor(all$effect_allele.exposure, levels=lev2)

# Quitamos espacios en blanco de la variable de alelo de exposición
all$effect_allele.exposure <- gsub(" ", "", all$effect_allele.exposure, fixed = TRUE)

# Multiplicamos por -1 todos los valores de beta de la exposición cuando el valor 
# el alelo de efecto de la exposición y del evento son diferentes.   
all$beta.exposure[all$effect_allele.exposure != all$effect_allele.outcome]  <- 
  all$beta.exposure[all$effect_allele.exposure != all$effect_allele.outcome] * -1

all[3,]

# 07 - Wald test ----------------------------------------------------------

# Corremos el test de Wald para efectos fijos y generamos un forest plot Igual
# que con la Aleatorización Mendeliana de varianza inversa ponderada (IVW) de
# efectos fijos

x <- all$beta.exposure                   # betas para los SNPs de la exposicion
sigmax <- all$se.exposure                # errores estandar de x
y <- all$beta.outcome                    # betas para los SNPs del evento
sigmay <- all$se.outcome                 # errores eestandar de y
all$Wald <- y / x                        # Wald estimate
all$Waldvar <- (sigmay^2 / x^2)          # using Burgess's method

dmres <- rma.uni(yi = all$Wald,          # Tamaño del efecto del evento
                 vi = all$Waldvar,       # Varianza de la muestra
                 slab = all$SNP,         # Etiquetas de los SNPs 
                 method = "FE"           # Modelo de efectos fijos
                 )
dmres

## 7.1
# Forest plot
forest(dmres, 
       atransf = exp,                    # transformación del los valores en el eje x, Exponencial
       xlab = " ",                       # Quitamos la etiqueta del eje x
       mlab = "Coronary Artery disease (OR)",  # Etiqueta de la columnna de SNPs 
       at = log(c(.5, 1,2)),             # Valores del eje x   
       xlim = c(-1.7,1.3),               # Rango del eje x
       cex = .8                          # Factor de expansion de los contenidos gel plot
       )

# Especificamos dimensiones de plot y guardamos el archivo 
par(mfrow = c(2,2)) 
par(mar = c(4,2,0,2) + .1)
png(file = "Output/forestplot_sol_urate_CAD.png") # Se guarda en la carpeta Output
forest(dmres, 
       atransf = exp,
       xlab = " ", 
       mlab = "Coronary Artery disease (OR)", 
       at = log(c(.5, 1,2)),
       xlim = c(-1.7,1.3),
       cex = .8)
title("MR uric acid to CAD")
dev.off()

# Data save ---------------------------------------------------------------

# Guardamos la base de dato con la que realizaremos nuestros análisis 

# write_csv(all, "Output/Data/urate_CAD.csv")

# En caso de que no haya podido realizar el trabajo, se puede cargar la base de 
# datos para realizar los análisis

# all <- read.csv("Output/Data/urate_CAD.csv")

# 08 - MR Estimates -------------------------------------------------------

# NOTA: Se utiliza el paquete MendelianRandomization para este análisis. Varias
# de las funciones en este paquete comparten nombre con el paquete TwoSampleMR

## 8.1
# Obtenemos las estimaciones para enfermedad cardiaca isquémica (IHD) usando el
# paquete de aleatorizacion mendeliana con el modelo de efectos fijos.

# Imputacion y formateo de datos para su uso en la estimación causal
MRInputObject <-
  mr_input(all$beta.exposure, 
           all$se.exposure, 
           all$beta.outcome, 
           all$se.outcome)

# Inverse-variance weighted method
mr_ivw(MRInputObject, 
       model = "fixed")     # Modelo de Efectos fijos

## 8.2
# Obtenemos las estimaciones para enfermedad cardiaca isquémica usando el paquete
# de aleatorizacion mendeliana con el modelo de efectos aleatorios
mr_ivw(MRInputObject)

## 8.3
# Obtenemos la mediana ponderada y la estimación MR-Egger
mr_median(MRInputObject)
mr_egger(MRInputObject)


# 09 - Mendelian Randomization Pleiotropy ---------------------------------

# Corremos la prueba Mendelian Randomization Pleiotropy RESidual Sum and Outlier
# (MR-Presso), este método ayuda a identificar SNPs que componen la IHC que son
# outlayers

## 9.1
# Para eso, tenemos que cargar datos simulados incluidos en el paquete MRPRESSO
data("SummaryStats") # No es necesario asignarlos a un objeto 

# Método global MR-Presso
mr_presso(BetaOutcome = "Y_effect",           # Beta del evento 
          BetaExposure = "E1_effect",         # Beta de la exposicion 
          SdOutcome = "Y_se",                 # Error estandar del evento  
          SdExposure = "E1_se",               # Error estandar de la exposicion 
          OUTLIERtest = TRUE,                 # Realiza la prueba de outliers
          DISTORTIONtest = TRUE,              # Realiza la prueba de distorsión de la estimación causal 
          data = SummaryStats,                # Base de datos a usar 
          NbDistribution = 1000,              # Elementos a simular para calcular el valor P empírico
          SignifThreshold = 0.05)             # Umbral del valor de significancia  

## 9.2
# Llevamos a cabo la prueba MR-Presso en la base de datos de IHC
mr_presso(BetaOutcome = "beta.outcome",       # Beta del evento 
          BetaExposure = "beta.exposure",     # Beta de la exposicion 
          SdOutcome = "se.outcome",           # Error estandar del evento  
          SdExposure = "se.exposure",         # Error estandar de la exposicion 
          OUTLIERtest = TRUE,                 # Realiza la prueba de outliers
          DISTORTIONtest = TRUE,              # Realiza la prueba de distorsión de la estimación causal 
          data = all,                         # Base de datos generada en el workshop 
          NbDistribution = 1000,              # Elementos a simular para calcular el valor P empírico
          SignifThreshold = 0.05)             # Umbral del valor de significancia  


## 9.3
# Volvemos a generar el forest plot con el modelos de efectos aleatorios 
# En caso de que se den, podemos eliminar outlayers

x <- all$beta.exposure            # betas para los SNPs de la exposicion
sigmax <- all$se.exposure         # errores estandar de x
y <- all$beta.outcome             # betas para los SNPs del evento
sigmay <- all$se.outcome          # errores eestandar de y
all$Wald2 <- y/x                  # Wald estimate
all$Waldvar2 <- (sigmay^2/x^2)    # using Burgess's method
all$lab <- paste(all$SNP, all$gene, sep=" ")

# Podemos eliminar de la base de datos todos los SNPs que sean outlayers.
# En este caso no tenemos outlayers 
# all <- all[all$SNP != "rs190212799",] 

dmres <- rma.uni(yi = all$Wald,     # Tamaño del efecto del evento
                 vi = all$Waldvar,  # Varianza de la muestra 
                 slab = all$lab,    # Etiquetas de los SNPs 
                 method = "REML")   # Modelo de efectos aleatorios (restricted maximum likelihood estimator) 
dmres

forest(dmres, 
       atransf = exp,
       xlab = " ", 
       mlab = "Coronary Artery disease (OR)", 
       at = log(c(.5, 1,2)),
       xlim = c(-1.7,1.3),
       cex = .8)


# Copyright ---------------------------------------------------------------

# YEAR: 2023
# Copyright (c) 2023: 
#   Created by: Mariaelisa Graff, PhD
#   Translation by: Carlos González-Carballo, MsC 

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this Documentation for the Mendelian Randomization Workshop, to deal in the
# Documentation without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Documentation, and to permit persons to whom the Documentation is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Documentation.
#
# THE DOCUMENTTION IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# DOCUMENTATION.