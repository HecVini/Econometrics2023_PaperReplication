# This function loads the data to estimate the model in Dennis, Quintero and Sieg (2019) with 1 metro area.

a1loadinitialdata_1ma <- function(WPmetro) {
  YIk <- WPmetro$YIk
  VIk <- WPmetro$VIk
  PIk <- WPmetro$PIk
  Pop <- WPmetro$Pop
  Corr <- WPmetro$Corr
  
  return(list(YIk = YIk, VIk = VIk, PIk = PIk, Pop = Pop, Corr = Corr))
}
