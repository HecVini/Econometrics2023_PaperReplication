### 1. Load Packages
library(tidyverse) # Easily Installand Load the 'Tidyverse' 
library(lubridate) # Make Dealing with Dates a Little Easier 
library(janitor) # Simple Tools for Examining and Cleaning Dirty Data 
library(data.table) # Extension of `data.frame`
library(tidylog) # Logging for 'dplyr' and 'tidyr' Functions 
library(rjson) # JSON for R
library(openxlsx) #Read, Write and Edit xlsx Files
library(haven) #Open .xpt and Stata (.do) files
library(nleqslv) #Solve system of non linear equations
library(rootSolve) #Solve system of non equations
library(knitr) #Make tables in RMarkdown
library(rio) #Deal with every kind of data type 
library(R.matlab) #Open .mat files
library(fmtr) #Get dataset descriptions
library(stargazer) #Make summary table in .Rmd files
library(kableExtra) #Make tables in markdown

## 2. Open Data
setwd("C:/Users/vinic/OneDrive/Mestrado/4Tri/Econometrics/PaperPersonalFolder/DataPrograms/Data") #Data programs folder in my PC

# 2.1 Correlation between Rent and Income
corrrenty_datasets = c('corrrentyChi99_k.dta', 'corrrentyChi03_k.dta', 'corrrentyNYC99_k.dta', 'corrrentyNYC03_k.dta') #List with the datasets used to calculate the correlation between income and rent.
corrrentyChi99_k.R = import(corrrenty_datasets[1]) #Open the fist on the list
corrrentyChi03_k.R = import(corrrenty_datasets[2]) #Open the second on the list
corrrentyNYC99_k.R = import(corrrenty_datasets[3]) #Open the third on the list
corrrentyNYC03_k.R = import(corrrenty_datasets[4]) #Open the fourth on the list

# 2.2 Equity in Chicago
AHS_full_datasets = c('hud1999.xpt','hud2003.xpt') #AHS long microdata for 1999 and 2003. HUD = Departament of Housing and Urban Development, US government responsable for doing the survey
hud1999.R = import(AHS_full_datasets[1]) #Open the fist on the list. That is a huge dataset.
hud2003.R = import(AHS_full_datasets[2]) #Open the second on the list. That is a huge dataset.

chicago_equity_datasets = c('Chicago_equity_1999.dta','Chicago_equity_2003.dta') #AHS survey filtered for Chicago and its equity data
Chicago_equity_1999.R = import(chicago_equity_datasets[1]) #Open the fist on the list. A short dataset with the results from 2.D (check the README file)
Chicago_equity_2003.R = import(chicago_equity_datasets[2]) #Open the second

Chicago_ucequity_income_iter_datasets = c('Chicago_ucequity_income_iter_1999.dta','Chicago_ucequity_income_iter_2003.dta') #Datasets used on 2.D. They are huge because have the resutls of multiple iteratiosn
Chicago_ucequity_income_iter_1999.R = import(Chicago_ucequity_income_iter_datasets[1]) #Open the first dataset
Chicago_ucequity_income_iter_2003.R = import(Chicago_ucequity_income_iter_datasets[2]) #Open the second one

# 2.3 ACS 
ACSReg.R = import('ACSReg.dta') #ACS = American Community Survey. Broader survey that collects a wide range of demographic and socioeconomic data about the U.S. population

## 3. Explore Data
# 3.1 Correlation between Rent and Income

corrrentyChi99_k = corrrentyChi99_k.R #Create a different object to be munipulated
ncol_corrrentyChi99_k = dim(corrrentyChi99_k)[2] #Number of columns in this dataset
corrrentyChi99_k_colnames = colnames(corrrentyChi99_k) #Get the colnames from it in a vector

corrrentyChi99_k_descriptions = numeric(ncol_corrrentyChi99_k) #Create a blank vector with dimension = number of cols
for (j in c(1:ncol_corrrentyChi99_k)) {
  label = attributes(corrrentyChi99_k[[j]])
  label = label[[1]]
  corrrentyChi99_k_descriptions[j] = label
} #Loop to replace zero in description vector by the dataset description

columntypes = sapply(corrrentyChi99_k, class) #Get the coltypes of the whole dataset
corrrentyChi99_k_coltypes = numeric(ncol_corrrentyChi99_k) #Blank vector to be fullied
for (k in c(1:ncol_corrrentyChi99_k)) {
  type = columntypes[[k]]
  corrrentyChi99_k_coltypes[k] = type
} #Loop to replace zero in t coltype by each column type

corrrentyChi99_k_metadata = tibble(colnames = corrrentyChi99_k_colnames,
                                   coltypes = corrrentyChi99_k_coltypes,
                                   short_description = corrrentyChi99_k_descriptions) %>% clean_names() #create a dataframe to explain metadata

correntyChi99_k_long_descriptions = tibble(colnames = corrrentyChi99_k_colnames, long_description = NA)
correntyChi99_k_long_descriptions[1,2] = "1 = Owned, 2 = Rent"
correntyChi99_k_long_descriptions[2,2] = "Identifcador unico de cada casa da pesquisa. Para Chicago 1999, sao 1909 control numbers unicos"
correntyChi99_k_long_descriptions[3,2] = "Codigo da regiao metropolitana do censo americano, usando a divisao de 1980. Para Chicago, 1600"
correntyChi99_k_long_descriptions[4,2] = "Valor do aluguel (rent) das obsrevacoes dos que nao tem ownership. Ver com calma, tem valores repetidos no topo"
correntyChi99_k_long_descriptions[5,2] = "Valor de mercado (value) das observacoes de quem tem ownership. Ver com calma tambem. NA para quem aluga"
correntyChi99_k_long_descriptions[6,2] = "Renda de todos os membros da casa. Imagino que em dolares correntes. Notar que ha obs zeradas"
correntyChi99_k_long_descriptions[7,2] = "01: House, apartment, flat, 02: Mobile home with no permanent room added, 03: Mobile home with permanent room added, 04: Housing unit in non transient hotel, motel etc., 05: Housing unit in transient hotel, motel etc., 06: Housing unit rooming house, 07: Boat or recreational vehicle, 09: Other housing units "
correntyChi99_k_long_descriptions[8,2] = "Renda que limita cada income bin (quantil de renda). Nesse caso, sao 16 valores diferentes"
correntyChi99_k_long_descriptions[9,2] = "Valor da casa que limita cada value bin (quantil do valor da casa). Nesse caso, sao 16 valores diferentes tambem"
correntyChi99_k_long_descriptions[10,2] = "Aluguel medio e corrigido pela inflacao para quantis de renda (income bins) para quem aluga"
correntyChi99_k_long_descriptions[11,2] = "Aluguel medio e corrigido pela inflacao para quantis de renda (income bins) por tipo para k-means (clusters)"
correntyChi99_k_long_descriptions[12,2] = "Valor medio das casas e corrigido pela inflacao para quantis de renda (income bins) para quem tem ownership"
correntyChi99_k_long_descriptions[13,2] = "Valor medio das casas e corrigido pela inflacao para quantis de renda (income bins) por tipo para k-means (clusters)"
correntyChi99_k_long_descriptions[14,2] = "Media da coluna 10 (meantentf), que e sao os mesmos valores (acho que existe pra substituir onde tem NA)"
correntyChi99_k_long_descriptions[15,2] = "UC = User Cost. Tirei da fig 1 do paper. Sao 16 valores unicos. Descobrir com calma o que e"
correntyChi99_k_long_descriptions[16,2] = "Media da coluna 11, meanrentfk. Pelo o que entendi, serve para colocar o valor medio do cluster para as observacoes NA"
correntyChi99_k_long_descriptions[17,2] = "UC para os clusters. Sao 48 valores unicos"
correntyChi99_k_long_descriptions[18,2] = "Ano. 1999 ou 2003"
correntyChi99_k_long_descriptions[19,2] = "Cidade. NYC ou Chicago"

corrrentyChi99_k_metadata = full_join(correntyChi99_k_long_descriptions,corrrentyChi99_k_metadata,by = "colnames",keep = FALSE) #Properly merge both dataframes
corrrentyChi99_k_metadata = corrrentyChi99_k_metadata %>% subset(select = c(colnames,coltypes,short_description,long_description)) #Remove duplicate column (make it in a better way)

# 3.2 Equity in Chicago
#hud1999 = hud1999.R #Create a different object to be manipulated
#hud2003 = hud2003.R 
AHS_codebook_link = 'https://www.census.gov/data-tools/demo/codebook/ahs/ahsdict.html?s_keyword=&s_variablelist='

Chicago_equity_1999 = Chicago_equity_1999.R 
Chicago_equity_2003 = Chicago_equity_2003.R 

Chicago_equity_1999 = Chicago_equity_1999.R #Create a different object to be munipulated
ncol_Chicago_equity_1999 = dim(Chicago_equity_1999)[2] #Number of columns in this dataset
Chicago_equity_1999_colnames = colnames(Chicago_equity_1999) #Get the colnames from it in a vector

columntypes = sapply(Chicago_equity_1999, class) #Get the coltypes of the whole dataset
Chicago_equity_1999_coltypes = numeric(ncol_Chicago_equity_1999) #Blank vector to be fullied
for (k in c(1:ncol_Chicago_equity_1999)) {
  type = columntypes[[k]]
  Chicago_equity_1999_coltypes[k] = type
} #Loop to replace zero in t coltype by each column type

Chicago_equity_1999_metadata = tibble(colnames = Chicago_equity_1999_colnames,
                                      coltypes = Chicago_equity_1999_coltypes) %>% clean_names() #create a dataframe to explain metadata
Chicago_equity_1999_metadata = Chicago_equity_1999_metadata %>%
  mutate(long_description = case_when(
    colnames == 'control' ~ "Control number, unique ID for each observation of the sample",
    colnames == 'otpinr' ~ "Outstanding Principal, valor do financiamento menos o que ja foi pago",
    colnames == 'equity' ~ "Valor de mercado da casa menos a hipoteca a ser paga (outstanding principal). Note que equity + otpinr = value" ))

Chicago_ucequity_income_iter_1999 = Chicago_ucequity_income_iter_1999.R

# 3.3 ACS
ACSReg = ACSReg.R
dim(ACSReg)
ncol_ACS = dim(ACSReg)[2] #Number of columns in this dataset
ACS_colnames = colnames(ACSReg) #Get the colnames from it in a vector

ACS_descriptions = numeric(ncol_ACS) #Create a blank vector with dimension = number of cols
for (j in c(1:ncol_ACS)) {
  label = attributes(ACSReg[[j]])
  label = label[[1]]
  ACS_descriptions[j] = label
} #Loop to replace zero in description vector by the dataset description

columntypes = sapply(ACSReg, class) #Get the coltypes of the whole dataset
ACS_coltypes = numeric(ncol_ACS) #Blank vector to be fullied
for (k in c(1:ncol_ACS)) {
  type = columntypes[[k]]
  ACS_coltypes[k] = type
} #Loop to replace zero in t coltype by each column type

ACS_metadata = tibble(colnames = ACS_colnames,
                      coltypes = ACS_coltypes,
                      short_description = ACS_descriptions) %>% clean_names() #create a dataframe to explain metadata
ACS_codebook_link = 'https://www2.census.gov/programs-surveys/acs/tech_docs/pums/data_dict/PUMS_Data_Dictionary_2022.pdf'

# 3.4 Matlab files
MetroAreasData.R = R.matlab::readMat('MetroAreasData.mat') #Read .mat file
MetroAreasData = MetroAreasData.R
MetroAreasData1 = MetroAreasData[[1]] #First subfolder of the file
MetroAreasData2 = MetroAreasData[[2]]
MetroAreasData3 = MetroAreasData[[3]]
MetroAreasData4 = MetroAreasData[[4]]
MetroAreasData5 = MetroAreasData[[5]]
MetroAreasData6 = MetroAreasData[[6]]
