# Article-Thrips---Austral-Enotmology
Script
library(rsq)
library(lme4)
library(vegan)
library(fitdistrplus)
library(dplyr)
library(glmmTMB)
library(car)
library(MASS)
library(interactions)
library(DHARMa)
library(MuMIn)
library(effectsize)
library(ggplot2)
library(jtools)
library(performance)
library(see)
library(qqplotr)
library(blmeco)
library(sjPlot)
library(bbmle) 
library(optimx)
library(nloptr)
library(dfoptim)
library(tidyr)
library(emmeans)
library(magrittr)
library(ggeffects)
library(viridis)
library(ggforce)
library(ggdist)
library(gghalves)
library(cowplot)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(broom.mixed)
library(jtools)
library(stargazer)
library(rsq)
library(lme4)
library(vegan)
library(fitdistrplus)
library(dplyr)
library(car)
library(jtools)
library(interactions)
library(performance)
library(see)
library(qqplotr)
library(blmeco)
library(sjPlot)
library(effects)
library(ggeffects)
setwd("C:/Users/Zoo/Documents/Daniel Ants")
caminho_arquivo <- "Abundancia_insetos_associados.csv"
Abundancia_insetos_associados <- read.csv(caminho_arquivo, header = T,  sep=";")

Abundancia_insetos_associados$blocos<-as.factor(Abundancia_insetos_associados$blocos)
Abundancia_insetos_associados$Treatment<-as.factor(Abundancia_insetos_associados$Treatment)
levels(Abundancia_insetos_associados$Treatment)
levels(Abundancia_insetos_associados$blocos)


#Modelos:
#Haplothrips
#com inflores em peso:
Haplothripes.sp. <- glmer(
  Haplothripes.sp. ~ Treatment + (1 | blocos),  
  family = poisson,
  data = Abundancia_insetos_associados,
  weights = Num_inflores.,
  glmerControl(optimizer = c("Nelder_Mead"))
)

residuos1 <- simulateResiduals(fittedModel = Haplothripes.sp., n=1000)
plot(residuos1)
summary(Haplothripes.sp.)
Anova(Haplothripes.sp.)

r.squaredGLMM(Haplothripes.sp.)
round(r.squaredGLMM(Haplothripes.sp.)*100,2)

# Extraindo a variância dos efeitos aleatórios
random_effect_variance <- VarCorr(Haplothripes.sp.)
variance_value <- as.numeric(random_effect_variance$blocos[1,1])

# Exibindo o modelo com a variância dos efeitos aleatórios
stargazer(Haplothripes.sp., 
          type = "text",
          style = "all",
          ci = TRUE,
          ci.level = 0.95,
          single.row = TRUE,
          p.auto = FALSE,
          digits = 3,
          model.numbers = FALSE,
          report = "vcsp*",
          out = "Premisses.Haplothrips.adult.html")

# Agora, podemos adicionar a variância da variável aleatória manualmente no arquivo de saída ou printá-la
cat("\n\nVariância dos efeitos aleatórios (blocos):", variance_value, "\n", file = "Premisses.Haplothrips.adult.html", append = TRUE)

# Extração das previsões com intervalos de confiança
predicoes <- emmeans(Haplothripes.sp., ~ Treatment)
predicoes_df <- as.data.frame(confint(predicoes)) # Inclui os intervalos de confiança

# Visualize as primeiras linhas para conferir as colunas
print(predicoes_df)

# Gráfico das previsões ponderadas
ggplot(predicoes_df, aes(x = Treatment, y = emmean)) +
  geom_point(size = 4, color = "darkblue") + # Pontos das médias preditas
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "darkblue") + # Intervalos de confiança
  geom_text(aes(label = "*"), nudge_y = 0.5, size = 15, color = "red") + # Opcional, para indicar significância
  labs(
    title = "b",
    x = "Treatment",
    y = expression(paste("Predicted abundance of ", italic(" Haplothripes fiebrigi"), " adult"))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),  # Tamanho e estilo do título
    axis.title.x = element_text(size = 14),              # Tamanho do título do eixo x
    axis.title.y = element_text(size = 14),              # Tamanho do título do eixo y
    axis.text.x = element_text(size = 12),               # Tamanho do texto do eixo x
    axis.text.y = element_text(size = 12)                # Tamanho do texto do eixo y
  )


#imaturos
#com inflores em peso:
imaturos.Haplothripes <- glmer(
  imaturos.Haplothripes ~ Treatment + (1 | blocos),  
  family = poisson,
  data = Abundancia_insetos_associados,
  weights = Num_inflores.,
  glmerControl(optimizer = c("bobyqa"))
)

residuos1 <- simulateResiduals(fittedModel = imaturos.Haplothripes, n=1000)
plot(residuos1)
summary(imaturos.Haplothripes)
Anova(imaturos.Haplothripes)

r.squaredGLMM(imaturos.Haplothripes)
round(r.squaredGLMM(imaturos.Haplothripes)*100,2)

# Extraindo a variância dos efeitos aleatórios
random_effect_variance <- VarCorr(imaturos.Haplothripes)
variance_value <- as.numeric(random_effect_variance$blocos[1,1])

# Exibindo o modelo com a variância dos efeitos aleatórios
stargazer(imaturos.Haplothripes, 
          type = "text",
          style = "all",
          ci = TRUE,
          ci.level = 0.95,
          single.row = TRUE,
          p.auto = FALSE,
          digits = 3,
          model.numbers = FALSE,
          report = "vcsp*",
          out = "Premisses.html")

# Agora, podemos adicionar a variância da variável aleatória manualmente no arquivo de saída ou printá-la
cat("\n\nVariância dos efeitos aleatórios (blocos):", variance_value, "\n", file = "Premisses.html", append = TRUE)

# Extração das previsões com intervalos de confiança
predicoes <- emmeans(imaturos.Haplothripes, ~ Treatment)
predicoes_df <- as.data.frame(confint(predicoes)) # Inclui os intervalos de confiança

# Visualize as primeiras linhas para conferir as colunas
print(predicoes_df)

# Gráfico das previsões ponderadas
ggplot(predicoes_df, aes(x = Treatment, y = emmean)) +
  geom_point(size = 4, color = "darkblue") + # Pontos das médias preditas
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "darkblue") + # Intervalos de confiança
  geom_text(aes(label = "*"), nudge_y = 0.5, size = 15, color = "red") + # Opcional, para indicar significância
  labs(
    title = "b",
    x = "Treatment",
    y = expression(paste("Predicted abundance of ", italic(" Haplothripes fiebrigi"), " immature"))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),  # Tamanho e estilo do título
    axis.title.x = element_text(size = 14),              # Tamanho do título do eixo x
    axis.title.y = element_text(size = 14),              # Tamanho do título do eixo y
    axis.text.x = element_text(size = 12),               # Tamanho do texto do eixo x
    axis.text.y = element_text(size = 12)                # Tamanho do texto do eixo y
  )


#Frankniella
#com inflores em peso:
Frankiniella.sp. <- glmer(
  Frankiniella.sp. ~ Treatment + (1 | blocos),  
  family = poisson,
  data = Abundancia_insetos_associados,
  weights = Num_inflores.,
  glmerControl(optimizer = c("bobyqa"))
)
residuos1 <- simulateResiduals(fittedModel = Frankiniella.sp., n=1000)
plot(residuos1)
summary(Frankiniella.sp.)
Anova(Frankiniella.sp.)

r.squaredGLMM(Frankiniella.sp.)
round(r.squaredGLMM(Frankiniella.sp.)*100,2)

# Extraindo a variância dos efeitos aleatórios
random_effect_variance <- VarCorr(Frankiniella.sp.)
variance_value <- as.numeric(random_effect_variance$blocos[1,1])

# Exibindo o modelo com a variância dos efeitos aleatórios
stargazer(Frankiniella.sp., 
          type = "text",
          style = "all",
          ci = TRUE,
          ci.level = 0.95,
          single.row = TRUE,
          p.auto = FALSE,
          digits = 3,
          model.numbers = FALSE,
          report = "vcsp*",
          out = "Premisses.Frankiniella.sp.adult.html")

# Agora, podemos adicionar a variância da variável aleatória manualmente no arquivo de saída ou printá-la
cat("\n\nVariância dos efeitos aleatórios (blocos):", variance_value, "\n", file = "Premisses.Frankiniella.sp.adult.html", append = TRUE)

# Extração das previsões ponderadas com intervalos de confiança
predicoes <- emmeans(Frankiniella.sp., ~ Treatment)
predicoes_df <- as.data.frame(summary(predicoes)) # Inclui os intervalos de confiança corretamente

# Visualize as primeiras linhas para conferir as colunas
print(predicoes_df)

# Gráfico das previsões ponderadas
ggplot(predicoes_df, aes(x = Treatment, y = emmean)) +
  geom_point(size = 4, color = "darkblue") + # Pontos das médias preditas
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "darkblue") + # Intervalos de confiança
  geom_text(aes(label = "*"), nudge_y = 0.5, size = 15, color = "red") + # Opcional, para indicar significância
  labs(
    title = "c",
    x = "Treatment",
    y = expression(paste("Predicted abundance of ", italic(" Frankiniella gemina"), " adults"))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),  # Tamanho e estilo do título
    axis.title.x = element_text(size = 14),              # Tamanho do título do eixo x
    axis.title.y = element_text(size = 14),              # Tamanho do título do eixo y
    axis.text.x = element_text(size = 12),               # Tamanho do texto do eixo x
    axis.text.y = element_text(size = 12)                # Tamanho do texto do eixo y
  )

#imaturos frankniella
#com inflores em peso:
imaturos.Frankiniella <- glmer(
  imaturos.Frankiniella ~ Treatment + (1 | blocos),  
  family = poisson,
  data = Abundancia_insetos_associados,
  weights = Num_inflores.,
  glmerControl(optimizer = c("bobyqa"))
)

residuos1 <- simulateResiduals(fittedModel = imaturos.Frankiniella, n=1000)
plot(residuos1)
summary(imaturos.Frankiniella) #infleunciado pelo tratamento
Anova(imaturos.Frankiniella)

r.squaredGLMM(imaturos.Frankiniella)
round(r.squaredGLMM(imaturos.Frankiniella)*100,2)

# Extraindo a variância dos efeitos aleatórios
random_effect_variance <- VarCorr(imaturos.Frankiniella)
variance_value <- as.numeric(random_effect_variance$blocos[1,1])

# Exibindo o modelo com a variância dos efeitos aleatórios
stargazer(imaturos.Frankiniella, 
          type = "text",
          style = "all",
          ci = TRUE,
          ci.level = 0.95,
          single.row = TRUE,
          p.auto = FALSE,
          digits = 3,
          model.numbers = FALSE,
          report = "vcsp*",
          out = "Premisses.imaturos.Frankiniella.html")

# Agora, podemos adicionar a variância da variável aleatória manualmente no arquivo de saída ou printá-la
cat("\n\nVariância dos efeitos aleatórios (blocos):", variance_value, "\n", file = "Premisses.imaturos.Frankiniella.html", append = TRUE)

# Extração das previsões ponderadas com intervalos de confiança
predicoes <- emmeans(imaturos.Frankiniella, ~ Treatment)
predicoes_df <- as.data.frame(summary(predicoes)) # Inclui os intervalos de confiança corretamente

# Visualize as primeiras linhas para conferir as colunas
print(predicoes_df)

# Gráfico das previsões ponderadas
ggplot(predicoes_df, aes(x = Treatment, y = emmean)) +
  geom_point(size = 4, color = "darkblue") + # Pontos das médias preditas
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, color = "darkblue") + # Intervalos de confiança
  geom_text(aes(label = "*"), nudge_y = 0.5, size = 15, color = "red") + # Opcional, para indicar significância
  labs(
    title = "d",
    x = "Treatment",
    y = expression(paste("Predicted abundance of ", italic(" Frankiniella gemina"), " immature"))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),  # Tamanho e estilo do título
    axis.title.x = element_text(size = 14),              # Tamanho do título do eixo x
    axis.title.y = element_text(size = 14),              # Tamanho do título do eixo y
    axis.text.x = element_text(size = 12),               # Tamanho do texto do eixo x
    axis.text.y = element_text(size = 12)                # Tamanho do texto do eixo y
  )


###################################################
setwd("C:/Users/Zoo/Documents/Daniel Ants")
caminho_arquivo <- "AssociadosXants number.csv"
AssociadosXants.number <- read.csv(caminho_arquivo, header = T,  sep=";")

#COMO A ABUNDANCIA DE FORMIGAS INFLUENCIA SOBRE ISSO:
shapiro.test(AssociadosXants.number$imaturos.Frankiniella)
imaturos.Frankiniella <- glmer(imaturos.Frankiniella ~ number.ants+(1|blocos), 
                               family=poisson,data = AssociadosXants.number,
                               glmerControl(optimizer = c("Nelder_Mead")))
residuos1 <- simulateResiduals(fittedModel = imaturos.Frankiniella, n=1000)
plot(residuos1)
summary(imaturos.Frankiniella)
Anova(imaturos.Frankiniella)

shapiro.test(AssociadosXants.number$imaturos.Haplothripes)
imaturos.Haplothripes <- glmer(imaturos.Haplothripes ~ number.ants+(1|blocos), 
                               family=poisson,data = AssociadosXants.number,
                               glmerControl(optimizer = c("Nelder_Mead")))
residuos1 <- simulateResiduals(fittedModel = imaturos.Haplothripes, n=1000)
plot(residuos1)
summary(imaturos.Haplothripes)
Anova(imaturos.Haplothripes)


#PARA OS INDIVIDUOS ADULTOS:
shapiro.test(AssociadosXants.number$Haplothripes.sp.)
Haplothripes.sp. <- glmer.nb(Haplothripes.sp. ~ number.ants+(1|blocos), 
                              data = AssociadosXants.number,
                          glmerControl(optimizer = c("bobyqa"), optCtrl = list(maxfun = 200000000)))
residuos1 <- simulateResiduals(fittedModel = Haplothripes.sp., n=1000)
plot(residuos1) #PROBLEMAS
summary(imaturos.Haplothripes)
Anova(imaturos.Haplothripes)

shapiro.test(AssociadosXants.number$Frankiniella.sp.)
Frankiniella.sp. <- glmer(Frankiniella.sp. ~ number.ants+(1|blocos), 
                          family=poisson,data = AssociadosXants.number,
                          glmerControl(optimizer = c("Nelder_Mead")))
residuos1 <- simulateResiduals(fittedModel = Frankiniella.sp., n=1000)
plot(residuos1)
summary(Frankiniella.sp.)
Anova(Frankiniella.sp.)

stargazer(Frankiniella.sp., Haplothripes.sp., type = "text",style= "all",ci=T,
          ci.level=0.95, single.row=TRUE, p.auto=F, digits=3,
          model.numbers=F,report= "vcsp*", out = "Premisses.html")
stargazer(imaturos.Haplothripes, imaturos.Frankiniella, type = "text",style= "all",ci=T,
          ci.level=0.95, single.row=TRUE, p.auto=F, digits=3,
          model.numbers=F,report= "vcsp*", out = "Premisses.html")

#Resumo; OS IMATUROS SÃO INFLUENCIADOS PELA PRESENÃA OU AUSENSIA DE FORMIGAS.
#NO CASO DOS Haplothrips immaturos, quando se tem formigas a abundancia Ã© maior
#Isso Ã© o oposto para Franklinella que diminui a  abundancia quando tem formigas.

#Entretanto esses dois resultados anteirores sÃ£o indeopendentes da abundancia de formigas,
#ou seja, se tem formigas na planta, eles diminuem em abundancia.

#Para os adultos nÃ£o encontramos efeitos dos tratamentos, entÃ£o se tem formigas,
#eles estarÃ£o lÃ¡. Porem, no caso de Haplothrips, ela prefere plantas com pouca inflorescencia.
