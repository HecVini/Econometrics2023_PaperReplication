# This function loads the data to estimate the model in Dennis, Quintero and Sieg (2019) with 1 metro area. 

a3loadinitialdata_2ma <- function(WPmetro, WPmetrom, WPAgr) {
  YIk <- WPmetro$YIk
  VIk <- WPmetro$VIk
  PIk <- WPmetro$PIk
  Pop <- WPmetro$Pop
  
  YIkm <- WPmetrom$YIk
  VIkm <- WPmetrom$Vik
  PIkm <- WPmetrom$PIk
  Popm <- WPmetrom$Pop
  
  Corra <- WPAgr$Corr
  
  return(list(YIk = YIk, VIk = VIk, PIk = PIk, Pop = Pop, 
              YIkm = YIkm, VIkm = VIkm, PIkm = PIkm, Popm = Popm, Corra = Corra))
}
