calculate.exclusivity <- function(df) {
  df.A <- df %>% dplyr::filter(group == "A_enriched")
  df.B <- df %>% dplyr::filter(group == "B_enriched")
  df.C <- df %>% dplyr::filter(group == "C_enriched")
  df.D <- df %>% dplyr::filter(group == "D_enriched")
  df.rest <- df %>% dplyr::filter(group == "rest_enriched")
  
  value.A <- (dim(dplyr::filter(df.A, exclusivity != "non-exclusive"))[1]/dim(df.A)[1])*100
  value.B <- (dim(dplyr::filter(df.B, exclusivity != "non-exclusive"))[1]/dim(df.B)[1])*100
  value.C <- (dim(dplyr::filter(df.C, exclusivity != "non-exclusive"))[1]/dim(df.C)[1])*100
  value.D <- (dim(dplyr::filter(df.D, exclusivity != "non-exclusive"))[1]/dim(df.D)[1])*100
  value.rest <- (dim(dplyr::filter(df.rest, exclusivity != "non-exclusive"))[1]/dim(df.rest)[1])*100
  
  df.out <- data.frame(percentage = c(value.A, value.B, value.C, value.D, value.rest), group = c("A", "B", "C", "D", "R"))
  
  return(df.out)
}
