#Estima um modelo apenas com os dados de Chi e estima outro com os dados de ambas as regiões e suas correlações

#setwd('')
WP<-readRDS('MetroAreasData.R')

a2<-a2Estimate(WP$Chi$YIk,WP$Chi$VIk,WP$Chi$PIk,WP$Chi$Pop,WP$Chi$Corr)

a3<-a3Estimate(WP$Chi$YIk,WP$Chi$VIk,WP$Chi$PIk,WP$Chi$Pop,WP$NYC$YIk,WP$NYC$VIk,WP$NYC$PIk,WP$NYC$Pop,WP$Agr$Corr)
