library(gglorenz)
setwd("C:/Users/vinic/OneDrive/Mestrado/4Tri/Econometrics/PaperPersonalFolder/Charts") #Create a chart subfolder to save all charts

# 4. Charts
#4.1 Income distribution in Chicago in 1999
Chicago99_income_distribution = corrrentyChi99_k %>% subset(select = c(control,zinc2f)) %>% tibble() #Subset control numbers and income data for each obsrvation
Chicago99_income_distribution_chart = 
  ggplot(data = Chicago99_income_distribution, aes(x = zinc2f)) +
  geom_histogram(fill = '#6c757d', color = "#000000",bins = 50) +
  scale_y_continuous(breaks = seq(0,360,40),labels = c('',seq(40,360,40))) +
  scale_x_continuous(breaks = seq(0,550000,50000),limits = c(0,550000),labels = seq(0,550,50)) +
  geom_hline(yintercept = 0, color = '#000000', size = 1) +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(size = 10, angle = 0, hjust = .5, vjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "#adb5bd", linewidth = 0.2),
        axis.ticks.y = element_blank())
ggsave('Chicago99_income_distribution.png', Chicago99_income_distribution_chart, 
       device = 'png', scale = 1, width = 15, units = "cm", dpi = 300, limitsize = TRUE)

#4.1.2 Income + Equity in Chicago in 1999 (ignorar)
Chicago99_income = corrrentyChi99_k %>% subset(select = c(control,zinc2f)) #Subset of incomes indexed by control
Chicago99_equity = Chicago_equity_1999 %>% subset(select = c(control,otpinr,equity)) #Subset of equity indexed by control
Chicago99_income_equity = full_join(Chicago99_income,Chicago99_equity, by = 'control') #Merge by control column
Chicago99_income_equity = Chicago99_income_equity %>% drop_na() #Remove NAS
Chicago99_income_equity = Chicago99_income_equity %>% tibble() #Make it a tibble
Chicago99_income_equity = Chicago99_income_equity %>% mutate(income_plus_equity = zinc2f + equity) #Create a column for the sum of income and equity 
Chicago99_income_equity = Chicago99_income_equity %>% filter(income_plus_equity >= 0) #Remove negatives #Create the difference between income and equity
Chicago99_income_equity = Chicago99_income_equity %>% mutate(diff_income_equity = income_plus_equity - zinc2f)

Chicago99_income_equity = full_join(Chicago99_income_distribution,Chicago99_income_equity, by = 'control') #Merge datasets by control number
Chicago99_income_equity = Chicago99_income_equity %>% subset(select = c(control,zinc2f.x,equity,income_plus_equity,diff_income_equity)) #Clean the dataset. Just some manipulaition
Chicago99_income_equity = Chicago99_income_equity %>% setnames(c('control','zinc2f','equity','income_plus_equity','diff_income_equity')) #Rename columns properly
Chicago99_income_equity = Chicago99_income_equity %>% drop_na(equity,zinc2f) #Remove NAs
Chicago99_income_equity = Chicago99_income_equity %>% arrange(desc(income_plus_equity)) #Order it

# 4.2 Total Equity
Chicago99_equity = Chicago_equity_1999 %>% subset(select = c(control,otpinr,equity)) %>% mutate(year = 1999) #Subset of equity indexed by control
Chicago03_equity = Chicago_equity_2003 %>% subset(select = c(control,otpinr,equity)) %>% mutate(year = 2003) #Subset of equity indexed by control

Chicago99_income_equity_chart = 
  ggplot(Chicago99_income_equity, aes(x = zinc2f, y = equity)) +
  geom_point(size = 1, color = '#adb5bd', alpha = 1) +
  geom_smooth(span = .75,color = '#1d3557') +
  scale_x_continuous(breaks = seq(0,400000,50000),limits = c(0,400000), labels = seq(0,400,50)) +
  scale_y_continuous(breaks = seq(0,450000,50000),labels = c('',seq(50,450,50)),limits = c(0,450000)) +
  geom_hline(yintercept = 0, color = '#000000', size = 1) +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(size = 10, angle = 0, hjust = .5, vjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "#adb5bd", linewidth = 0.2),
        axis.ticks.y = element_blank()) #Chart of distribution of income vs equity
ggsave('Chicago99_income_equity.png', Chicago99_income_equity_chart, 
       device = 'png', scale = 1, width = 15, units = "cm", dpi = 300, limitsize = TRUE)

## 4.3 Income Distribution
# 4.3.1 Income Lorenz
Chicago99_income_deciles = corrrentyChi99_k %>% subset(select = c(control,zinc2f)) %>% filter(zinc2f > 0)%>% tibble() #Make a tibble to make the lorenz curve with nonzero income
Chicago99_income_deciles = Chicago99_income_deciles %>% 
  mutate(income_decile = ntile(zinc2f, 10)) #Categorize each obs by income decile
Chicago99_income_deciles = Chicago99_income_deciles %>% arrange(zinc2f) #Order by income

Chicago99_income_deciles_totals = Chicago99_income_deciles %>%
  group_by(income_decile) %>% summarise(decile_income = sum(zinc2f)) #See total income by decile
Chicago99_income_deciles_totals = Chicago99_income_deciles_totals %>%
  mutate(cumulative_income = cumsum(decile_income)) #Check cumulative income by decile
Chicago99_income_deciles_totals = Chicago99_income_deciles_totals %>% mutate(income_share = cumulative_income/sum(decile_income)) #Check cumulative income share of each decile. E.g: share of D3 is (D1 + D2 + D3)/all income

Chicago99_income_deciles_chart = 
  ggplot(Chicago99_income_deciles,aes(x = zinc2f)) +
  stat_lorenz(geom = "area", alpha = 0.9,desc = FALSE,fill = '#1d3557') +
  scale_x_continuous(breaks = seq(0,1,.1),labels = c('','P10','P20','P30','P40','P50','P60','P70','P80','P90','')) +
  scale_y_continuous(breaks = seq(0,1,.1),labels = c('','10%','20%','30%','40%','50%','60%','70%','80%','90%','100%')) +
  geom_hline(yintercept = 0, color = '#000000', size = 1) +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(size = 10, angle = 0, hjust = .5, vjust = 1),
        axis.text.y = element_text(size = 10,hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "#adb5bd", linewidth = 0.2),
        axis.ticks.y = element_blank())
Chicago99_income_deciles_chart
ggsave('Chicago99_income_deciles_chart.png', Chicago99_income_deciles_chart, 
       device = 'png', scale = 1, width = 15, units = "cm", dpi = 300, limitsize = TRUE)


## 4.4 Rent (para depois)
Chicago99_rent = corrrentyChi99_k %>% subset(select = c(control,rentf)) %>% tibble() #Only select control nubmers and rent
Chicago99_rent = Chicago99_rent %>% drop_na() #Remove NAs
Chicago99_rent = Chicago99_rent %>% arrange(desc(rentf)) #Order by rent pay
