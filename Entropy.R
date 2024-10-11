##### Libraries ####
library( plotly )
library( caret )

##### Functions ####

MinMax <- function( x, new_min = 0, new_max = 1 ){
  out <- ( ( x - min(x) ) / ( max(x) - min(x) ) ) * ( new_max - new_min ) + new_min
  return( out )
}

Shannon <- function( p, ad = 1 ){
  out <- -apply( p * log(p,2), 1, sum )
  out[is.na(out)] <- 0 
  return( out )
}

Renyi <- function( p, Q = 1.1 ){ # Q >= 0, Q != 1 
  out <- 1 / (1 - Q) * apply( p ^ Q, 1, function(x){ log(sum(x),2) } )
  return( out )
}

Tsallis <- function( p, Q ){ # Q != 1 
  out <- 1 / (Q - 1) * ( 1 - apply( p ^ Q, 1, function(x){ sum(x) } ) )
  return( out )
}

Sharma_Mittal <- function( p, R = 1.5, Q = 1.5 ){ # R != 1, R != Q, 
  out <- ( 1 / (1-R) ) * ( apply( p ^ Q, 1, sum ) ^ ( ( 1 - R ) / ( 1 - Q ) ) )
  return( out )
}

Sharma_Taneja <- function( p, A = 1, B = 1 ){ # A B >=0, A != B, 
  out <- ( 2 ^ ( 1 - A ) - 2 ^ ( 1 - B ) ) ^ (-1) * ( apply( p ^ A, 1, sum ) - apply( p ^ B, 1, sum ) )
  return( out )
}

Kapur <- function( p, A = 1.5, B = 1 ){ # A != 1, A B >= 0, A + B -1 > 0, 
  out <- ( 1 / ( 1 - A ) ) * log( apply( p ^ ( A + B - 1 ), 1, sum ) / apply( p ^ B, 1, sum ) )
  return( out )
}

Gini <- function( p, sc = 2 ){
  out <- ( 1 - apply( p^2, 1, sum ) ) * sc # * 2 aby skala byla ta sama
  return( out )
}

Miss <- function( p, power = 1, sc = 2 ){
  out <- ( 1 - apply( p, 1, max ) ^ power ) * sc # * 2 aby skala byla ta sama
  return( out )
}

SimVar <- function( n = 1000, distType ){

  name <- names( distType )
  Data <- sapply( name, function( typ, dist ){ 
    eval( parse( text = sprintf( "%s(%s,%s)", dist[[ typ ]][[ 1 ]], n, 
                                 paste0( dist[[ typ ]][[ 2 ]], collapse = "," ) ) ) )
    }, distType )

  if( length(distType) > 1 ){

    Data <- apply( Data, 2, MinMax, -4, 4 )
    linProb <- apply( Data, 1, sum )
    normlinProb <- MinMax( linProb, -4, 4 )

  }else{
    
    linProb <- apply( Data, 1, sum )
    normlinProb <- MinMax( linProb, -4, 4 )
    
  }

  
  prob <- binomial()$linkinv( as.double( linProb ) ) # prob <= 0.5 -> 0 else 1
  probnorm <- binomial()$linkinv( as.double( normlinProb ) )
  
  out <- list( Probs = prob, ProbsNorm = probnorm, Quantile = linProb, normQuantile = normlinProb, Data = Data )
  attr( out, "types") <- name
  
  return( out )
  
}

Plots <- function( dat, slope1 = TRUE, slope2 = TRUE  ){
  
  Dist <- paste0( attr( dat, "types" ), collapse = "+" )

  #### Data ####  
  p <- sort( dat$Probs )
  p <- data.frame( p = p, q = 1 - p )
  p <- data.frame( p = p$p, Shannon = Shannon( p ), 
                   Gini = Gini( p, 1 ), Gini2 = Gini( p, 2 ), Miss = Miss( p, 1, 1 ), Miss2 = Miss( p, 1, 2 ),
                   Renyi0 = Renyi( p, 0 ), Renyi0.5 = Renyi( p, 0.5 ), Renyi1.5 = Renyi( p, 1.5 ), Renyi2.0 = Renyi( p, 2.0 ),
                   Tsallis0 = Tsallis( p, 0 ), Tsallis0.5 = Tsallis( p, 0.5 ), Tsallis1.5 = Tsallis( p, 1.5 ), Tsallis2.0 = Tsallis( p, 2.0 )
                  )

  p_norm <- sort( dat$ProbsNorm )
  p_norm <- data.frame( p = p_norm, q = 1 - p_norm )
  p_norm <- data.frame( p = p_norm$p, Shannon = Shannon( p_norm ), 
                        Gini = Gini( p_norm, 1 ), Gini2 = Gini( p_norm, 2 ), Miss = Miss( p_norm, 1, 1 ), Miss2 = Miss( p_norm, 1, 2 ),
                        Renyi0 = Renyi( p_norm, 0 ), Renyi0.5 = Renyi( p_norm, 0.5 ), Renyi1.5 = Renyi( p_norm, 1.5 ), Renyi2.0 = Renyi( p_norm, 2.0 ),
                        Tsallis0 = Tsallis( p_norm, 0 ), Tsallis0.5 = Tsallis( p_norm, 0.5 ), Tsallis1.5 = Tsallis( p_norm, 1.5 ), Tsallis2.0 = Tsallis( p_norm, 2.0 ) 
                       )
  
  q <- dat$Quantile
  q_norm <- dat$normQuantile
  
  FIG <- Slope <- list()
  
  #### Plot 1 ####  
  fit <- density( q )
  fig <- plot_ly( x = q, type = "histogram", name = "Histogram" )
  fig <- add_trace( fig, x = fit$x, y = fit$y, type = "scatter", mode = "lines", fill = "tozeroy", yaxis = "y2", name = "Density" )
  fig <- layout( fig, title = paste0( Dist, " : variable distribution" ),
                 xaxis = list( title = "Values" ),
                 yaxis2 = list( overlaying = "y", side = "right" ) )
  FIG[[ 1 ]] <- fig
  
  #### Plot 2 ####
  fit <- density( q_norm )
  fig <- plot_ly( x = q_norm, type = "histogram", name = "Histogram" )
  fig <- add_trace( fig, x = fit$x, y = fit$y, type = "scatter", mode = "lines", fill = "tozeroy", yaxis = "y2", name = "Density" )
  fig <- layout( fig, title = paste0( Dist, " :  normalized variable distribution" ),
                 xaxis = list( title = "Values" ),
                 yaxis2 = list( overlaying = "y", side = "right" ) )
  FIG[[ 2 ]] <- fig
  
  #### Plot 3 ####
  fit <- density( p$p )
  fig <- plot_ly( x = p$p, type = "histogram", name = "Histogram" )
  fig <- add_trace( fig, x = fit$x, y = fit$y, type = "scatter", mode = "lines", fill = "tozeroy", yaxis = "y2", name = "Density" )
  fig <- layout( fig, title = paste0( Dist, " - probability distribution from the model" ),
                 xaxis = list( title = TeX("p_\\text{2}") ),
                 yaxis2 = list( overlaying = "y", side = "right" ) )
  fig <- config( fig, mathjax = 'cdn' )
  FIG[[ 3 ]] <- fig
  
  #### Plot 4 ####
  fit <- density( p_norm$p )
  fig <- plot_ly( x = p_norm$p, type = "histogram", name = "Histogram" )
  fig <- add_trace( fig, x = fit$x, y = fit$y, type = "scatter", mode = "lines", fill = "tozeroy", yaxis = "y2", name = "Density" )
  fig <- layout( fig, title = paste0( Dist, " - probability distribution from the model" ),
                 xaxis = list( title = TeX("p_\\text{2}") ),
                 yaxis2 = list( overlaying = "y", side = "right" ) )
  fig <- config( fig, mathjax = 'cdn' )
  FIG[[ 4 ]] <- fig
  
  #### Plot 5 ####
  fig <- plot_ly( x = p$p, y = p[,"Shannon"], type = 'scatter', mode = 'markers', name = "Shannon" )
  fig <- layout( fig, title = paste0( Dist, " distribution" ), 
                 xaxis = list( title = "Probability", range=c( 0, 1 ) ), 
                 yaxis = list( title = "Measure" ) )
  for( i in 3:6 ){
    colName <- colnames( p )[i]
    fig <- add_trace( fig, x = p$p, y = p[,colName], type = 'scatter', mode = 'markers', name = colName )
  }
  
  if( slope1 ){
    p_ <- p[ p$p >= 0.5 & p$p <= 0.6, ]
    for( i in c("Shannon","Gini","Gini2","Miss","Miss2") ) {
      
      reg <- lm( formula( paste0( i, "~ p" ) ), data = p_ )

      pred <- data.frame( x = p$p, y = predict( reg, data.frame( p = p$p ) ) )
      pred <- pred[ pred$y <= 1.05, ]
      
      par <- ifelse( is.na( coef(reg)[2] ), 0, coef(reg)[2] )
      fig <- add_trace( fig, x = pred$x, y = pred$y, mode = 'lines', 
                        name = paste0( i, ': ', round( par, 2 ) ) )
      
    }
    
  }
  FIG[[ 5 ]] <- fig
  
  #### Plot 6 ####
  fig <- plot_ly( x = p$p, y = p[,"Shannon"], type = 'scatter', mode = 'markers', name = "Shannon" )
  fig <- layout( fig, title = paste0( Dist, " distribution" ), 
                 xaxis = list( title = "Probability", range=c( 0, 1 ) ), 
                 yaxis = list( title = "Measure" ) )
  for( i in 7:10 ){
    colName <- colnames( p )[i]
    fig <- add_trace( fig, x = p$p, y = p[,colName], type = 'scatter', mode = 'markers', name = colName )
  }
  
  if( slope1 ){
    p_ <- p[ p$p >= 0.5 & p$p <= 0.6, ]
    for( i in 7:10 ) {
      
      colName <- colnames( p )[i]
      reg <- lm( formula( paste0( colName, "~ p" ) ), data = p_ )
      
      pred <- data.frame( x = p$p, y = predict( reg, data.frame( p = p$p ) ) )
      pred <- pred[ pred$y <= 1.05, ]
      
      par <- ifelse( is.na( coef(reg)[2] ), 0, coef(reg)[2] )
      fig <- add_trace( fig, x = pred$x, y = pred$y, mode = 'lines', 
                        name = paste0( colName, ': ', round( par, 2 ) ) )
      
    }
  }
  FIG[[ 6 ]] <- fig
  
  #### Plot 7 ####
  fig <- plot_ly( x = p$p, y = p[,"Shannon"], type = 'scatter', mode = 'markers', name = "Shannon" )
  fig <- layout( fig, title = paste0( Dist, " distribution" ), 
                 xaxis = list( title = "Probability", range=c( 0, 1 ) ), 
                 yaxis = list( title = "Measure" ) )
  for( i in 11:14 ){
    colName <- colnames( p )[i]
    fig <- add_trace( fig, x = p$p, y = p[,colName], type = 'scatter', mode = 'markers', name = colName )
  }
  
  if( slope1 ){
    p_ <- p[ p$p >= 0.5 & p$p <= 0.6, ]
    for( i in 11:14 ) {
      
      colName <- colnames( p )[i]
      reg <- lm( formula( paste0( colName, "~ p" ) ), data = p_ )
      
      pred <- data.frame( x = p$p, y = predict( reg, data.frame( p = p$p ) ) )
      pred <- pred[ pred$y <= 1.05, ]
      
      par <- ifelse( is.na( coef(reg)[2] ), 0, coef(reg)[2] )
      fig <- add_trace( fig, x = pred$x, y = pred$y, mode = 'lines', 
                        name = paste0( colName, ': ', round( par, 2 ) ) )
      
    }
  }
  FIG[[ 7 ]] <- fig
  
  #### Plot 8 ####
  RQ <- expand.grid( AR = seq( 0, 2, by = 0.5 ), BQ = seq( 0, 2, by = 0.5 ) )
  
  RQ_ <- RQ[ RQ$AR != 1 & RQ$AR != RQ$BQ, ] 
  p_ <- data.frame( p = c(), R = c(), Q = c(), val = c() )
  slope <- c()
  for( i in 1:nrow(RQ_) ){
    p__ <- data.frame( p = p$p, R = RQ_$AR[i], Q = RQ_$BQ[i], 
                       val = Sharma_Mittal( data.frame( p = p$p, q = 1 - p$p ), RQ_$AR[i], RQ_$BQ[i] ) )
    p_ <- rbind( p_, p__ )
    p___ <- p__[ p__$p >= 0.5 & p__$p <= 0.6, ]
    reg <- lm( val ~ p, data = p___ )
    slope <- c( slope, round( coef(reg)[2], 2 ) )
  }
  Slope$Sharma_Mittal <- cbind( RQ_, slope )
  Slope$Sharma_Mittal <- reshape( Slope$Sharma_Mittal, idvar = "AR", timevar = "BQ", direction = "wide" )
  
  p_$Q <- as.factor( p_$Q )
  
  fig <- plot_ly( x = p_$p, y = p_$R, z = p_$val, color = p_$Q )
  fig <- add_markers( fig )
  fig <- layout( fig, title = list( text = paste0( Dist, " distribution" ), y = 0.8 ),
                 scene = list(xaxis = list(title = 'Probability', range=c( 0, 1 ) ),
                             yaxis = list(title = 'R'),
                             zaxis = list(title = 'Sharma-Mittal'),
                             camera = list( eye = list( x = 1.8, y = -1.5, z = 0.5 ) )
                             ),
                 legend = list( x = 1, y = 0.1 ), 
                 margin = list( l = 0, r = 0, b = 0, t = 0, pad = 0 ),
                 annotations = list( x = 1, y = 0.25, text = "Q", xref = 'paper', yref = 'paper', showarrow = FALSE )
                )

  FIG[[ 8 ]] <- fig
  
  #### Plot 9 ####
  RQ_ <- RQ[ RQ$AR >= 0 & RQ$BQ >= 0 & RQ$AR != RQ$BQ, ] 
  p_ <- data.frame( p = c(), R = c(), Q = c(), val = c() )
  slope <- c()
  for( i in 1:nrow(RQ_) ){
    p__ <- data.frame( p = p$p, R = RQ_$AR[i], Q = RQ_$BQ[i], 
                       val = Sharma_Taneja( data.frame( p = p$p, q = 1 - p$p ), RQ_$AR[i], RQ_$BQ[i] ) )
    p_ <- rbind( p_, p__ )
    p___ <- p__[ p__$p >= 0.5 & p__$p <= 0.6, ]
    reg <- lm( val ~ p, data = p___ )
    slope <- c( slope, round( coef(reg)[2], 2 ) )
  }
  Slope$Sharma_Taneja <- cbind( RQ_, slope )
  Slope$Sharma_Taneja <- reshape( Slope$Sharma_Taneja, idvar = "AR", timevar = "BQ", direction = "wide" )
  
  p_$Q <- as.factor( p_$Q )
  
  fig <- plot_ly( x = p_$p, y = p_$R, z = p_$val, color = p_$Q )
  fig <- add_markers( fig )
  fig <- layout( fig, title = list( text = paste0( Dist, " distribution" ), y = 0.8 ), 
                 scene = list(xaxis = list(title = 'Probability', range=c( 0, 1 ) ),
                              yaxis = list(title = 'A'),
                              zaxis = list(title = 'Sharma-Taneja'),
                              camera = list( eye = list( x = 1.8, y = -1.5, z = 0.5 ) )
                             ),
                 legend = list( x = 1, y = 0.1 ), 
                 margin = list( l = 0, r = 0, b = 0, t = 0, pad = 0 ),
                 annotations = list( x = 1, y = 0.25, text = "B", xref = 'paper', yref = 'paper', showarrow = FALSE )
  )
  FIG[[ 9 ]] <- fig
  
  #### Plot 10 ####
  RQ_ <- RQ[ RQ$AR != 1 & RQ$AR >= 0 & RQ$BQ >= 0 & (RQ$AR + RQ$BQ -1 > 0 ), ] 
  p_ <- data.frame( p = c(), R = c(), Q = c(), val = c() )
  slope <- c()
  for( i in 1:nrow(RQ_) ){
    p__ <- data.frame( p = p$p, R = RQ_$AR[i], Q = RQ_$BQ[i], 
                       val = Kapur( data.frame( p = p$p, q = 1 - p$p ), RQ_$AR[i], RQ_$BQ[i] ) )
    p_ <- rbind( p_, p__ )
    p___ <- p__[ p__$p >= 0.5 & p__$p <= 0.6, ]
    reg <- lm( val ~ p, data = p___ )
    slope <- c( slope, round( coef(reg)[2], 2 ) )
  }
  Slope$Kapur <- cbind( RQ_, slope )
  Slope$Kapur <- reshape( Slope$Kapur, idvar = "AR", timevar = "BQ", direction = "wide" )
  
  p_$Q <- as.factor( p_$Q )
  
  fig <- plot_ly( x = p_$p, y = p_$R, z = p_$val, color = p_$Q )
  fig <- add_markers( fig )
  fig <- layout( fig, title = list( text = paste0( Dist, " distribution" ), y = 0.8 ), 
                 scene = list(xaxis = list(title = 'Probability', range=c( 0, 1 ) ),
                              yaxis = list(title = 'A'),
                              zaxis = list(title = 'Kapur' ),
                              camera = list( eye = list( x = 1.8, y = -1.5, z = 0.5 ) )
                              ),
                 legend = list( x = 1, y = 0.1 ), 
                 margin = list( l = 0, r = 0, b = 0, t = 0, pad = 0 ),
                 annotations = list( x = 1, y = 0.25, text = "B", xref = 'paper', yref = 'paper', showarrow = FALSE )
  )
  
  FIG[[ 10 ]] <- fig
  
  #### Plot 11 ####
  p <- p_norm
  fig <- plot_ly( x = p$p, y = p[,"Shannon"], type = 'scatter', mode = 'markers', name = "Shannon" )
  fig <- layout( fig, title = paste0( Dist, " normalized distribution" ), 
                 xaxis = list( title = "Probability", range=c( 0, 1 ) ), 
                 yaxis = list( title = "Measure" ) )
  for( i in 3:6 ){
    colName <- colnames( p )[i]
    fig <- add_trace( fig, x = p$p, y = p[,colName], type = 'scatter', mode = 'markers', name = colName )
  }
  
  if( slope2 ){
    p_ <- p[ p$p >= 0.5 & p$p <= 0.6, ]
    for( i in c("Shannon","Gini","Gini2","Miss","Miss2") ) {
      
      reg <- lm( formula( paste0( i, "~ p" ) ), data = p_ )

      pred <- data.frame( x = p$p, y = predict( reg, data.frame( p = p$p ) ) )
      pred <- pred[ pred$y <= 1.05, ]
      
      par <- ifelse( is.na( coef(reg)[2] ), 0, coef(reg)[2] )
      fig <- add_trace( fig, x = pred$x, y = pred$y, mode = 'lines', 
                        name = paste0( i, ': ', round( par, 2 ) ) )
      
    }
    
  }
  FIG[[ 11 ]] <- fig
  
  #### Plot 12 ####
  p <- p_norm
  fig <- plot_ly( x = p$p, y = p[,"Shannon"], type = 'scatter', mode = 'markers', name = "Shannon" )
  fig <- layout( fig, title = paste0( Dist, " normalized distribution" ), 
                 xaxis = list( title = "Probability", range=c( 0, 1 ) ), 
                 yaxis = list( title = "Measure" ) )
  for( i in 7:10 ){
    colName <- colnames( p )[i]
    fig <- add_trace( fig, x = p$p, y = p[,colName], type = 'scatter', mode = 'markers', name = colName )
  }
  
  if( slope2 ){
    p_ <- p[ p$p >= 0.5 & p$p <= 0.6, ]
    for( i in 7:10 ) {
      
      colName <- colnames( p )[i]
      reg <- lm( formula( paste0( colName, "~ p" ) ), data = p_ )
      
      pred <- data.frame( x = p$p, y = predict( reg, data.frame( p = p$p ) ) )
      pred <- pred[ pred$y <= 1.05, ]
      
      par <- ifelse( is.na( coef(reg)[2] ), 0, coef(reg)[2] )
      fig <- add_trace( fig, x = pred$x, y = pred$y, mode = 'lines', 
                        name = paste0( colName, ': ', round( par, 2 ) ) )
      
    }
  }
  FIG[[ 12 ]] <- fig
  
  #### Plot 13 ####
  p <- p_norm
  fig <- plot_ly( x = p$p, y = p[,"Shannon"], type = 'scatter', mode = 'markers', name = "Shannon" )
  fig <- layout( fig, title = paste0( Dist, " normalized distribution" ), 
                 xaxis = list( title = "Probability", range=c( 0, 1 ) ), 
                 yaxis = list( title = "Measure" ) )
  for( i in 11:14 ){
    colName <- colnames( p )[i]
    fig <- add_trace( fig, x = p$p, y = p[,colName], type = 'scatter', mode = 'markers', name = colName )
  }
  
  if( slope2 ){
    p_ <- p[ p$p >= 0.5 & p$p <= 0.6, ]
    for( i in 11:14 ) {
      
      colName <- colnames( p )[i]
      reg <- lm( formula( paste0( colName, "~ p" ) ), data = p_ )
      
      pred <- data.frame( x = p$p, y = predict( reg, data.frame( p = p$p ) ) )
      pred <- pred[ pred$y <= 1.05, ]
      
      par <- ifelse( is.na( coef(reg)[2] ), 0, coef(reg)[2] )
      fig <- add_trace( fig, x = pred$x, y = pred$y, mode = 'lines', 
                        name = paste0( colName, ': ', round( par, 2 ) ) )
      
    }
  }
  FIG[[ 13 ]] <- fig
  
  #### Plot 14 ####
  p <- p_norm
  RQ <- expand.grid( AR = seq( 0, 2, by = 0.5 ), BQ = seq( 0, 2, by = 0.5 ) )
  
  RQ_ <- RQ[ RQ$AR != 1 & RQ$AR != RQ$BQ, ] 
  p_ <- data.frame( p = c(), R = c(), Q = c(), val = c() )
  slope <- c()
  for( i in 1:nrow(RQ_) ){
    p__ <- data.frame( p = p$p, R = RQ_$AR[i], Q = RQ_$BQ[i], 
                       val = Sharma_Mittal( data.frame( p = p$p, q = 1 - p$p ), RQ_$AR[i], RQ_$BQ[i] ) )
    p_ <- rbind( p_, p__ )
    p___ <- p__[ p__$p >= 0.5 & p__$p <= 0.6, ]
    reg <- lm( val ~ p, data = p___ )
    slope <- c( slope, round( coef(reg)[2], 2 ) )
  }
  Slope$Sharma_Mittal_norm <- cbind( RQ_, slope )
  Slope$Sharma_Mittal_norm <- reshape( Slope$Sharma_Mittal_norm, idvar = "AR", timevar = "BQ", direction = "wide" )
  
  p_$Q <- as.factor( p_$Q )
  
  fig <- plot_ly( x = p_$p, y = p_$R, z = p_$val, color = p_$Q )
  fig <- add_markers( fig )
  fig <- layout( fig, title = list( text = paste0( Dist, " normalized distribution" ), y = 0.8 ),
                 scene = list(xaxis = list(title = 'Probability', range=c( 0, 1 ) ),
                              yaxis = list(title = 'R'),
                              zaxis = list(title = 'Sharma-Mittal'),
                              camera = list( eye = list( x = 1.8, y = -1.5, z = 0.5 ) )
                 ),
                 legend = list( x = 1, y = 0.1 ), 
                 margin = list( l = 0, r = 0, b = 0, t = 0, pad = 0 ),
                 annotations = list( x = 1, y = 0.25, text = "Q", xref = 'paper', yref = 'paper', showarrow = FALSE )
  )
  
  FIG[[ 14 ]] <- fig
  
  #### Plot 15 ####
  p <- p_norm
  RQ_ <- RQ[ RQ$AR >= 0 & RQ$BQ >= 0 & RQ$AR != RQ$BQ, ] 
  p_ <- data.frame( p = c(), R = c(), Q = c(), val = c() )
  slope <- c()
  for( i in 1:nrow(RQ_) ){
    p__ <- data.frame( p = p$p, R = RQ_$AR[i], Q = RQ_$BQ[i], 
                       val = Sharma_Taneja( data.frame( p = p$p, q = 1 - p$p ), RQ_$AR[i], RQ_$BQ[i] ) )
    p_ <- rbind( p_, p__ )
    p___ <- p__[ p__$p >= 0.5 & p__$p <= 0.6, ]
    reg <- lm( val ~ p, data = p___ )
    slope <- c( slope, round( coef(reg)[2], 2 ) )
  }
  Slope$Sharma_Taneja_norm <- cbind( RQ_, slope )
  Slope$Sharma_Taneja_norm <- reshape( Slope$Sharma_Taneja_norm, idvar = "AR", timevar = "BQ", direction = "wide" )
  
  p_$Q <- as.factor( p_$Q )
  
  fig <- plot_ly( x = p_$p, y = p_$R, z = p_$val, color = p_$Q )
  fig <- add_markers( fig )
  fig <- layout( fig, title = list( text = paste0( Dist, " normalized distribution" ), y = 0.8 ), 
                 scene = list(xaxis = list(title = 'Probability', range=c( 0, 1 ) ),
                              yaxis = list(title = 'A'),
                              zaxis = list(title = 'Sharma-Taneja'),
                              camera = list( eye = list( x = 1.8, y = -1.5, z = 0.5 ) )
                 ),
                 legend = list( x = 1, y = 0.1 ), 
                 margin = list( l = 0, r = 0, b = 0, t = 0, pad = 0 ),
                 annotations = list( x = 1, y = 0.25, text = "B", xref = 'paper', yref = 'paper', showarrow = FALSE )
  )
  FIG[[ 15 ]] <- fig
  
  #### Plot 16 ####
  p <- p_norm
  RQ_ <- RQ[ RQ$AR != 1 & RQ$AR >= 0 & RQ$BQ >= 0 & (RQ$AR + RQ$BQ -1 > 0 ), ] 
  p_ <- data.frame( p = c(), R = c(), Q = c(), val = c() )
  slope <- c()
  for( i in 1:nrow(RQ_) ){
    p__ <- data.frame( p = p$p, R = RQ_$AR[i], Q = RQ_$BQ[i], 
                       val = Kapur( data.frame( p = p$p, q = 1 - p$p ), RQ_$AR[i], RQ_$BQ[i] ) )
    p_ <- rbind( p_, p__ )
    p___ <- p__[ p__$p >= 0.5 & p__$p <= 0.6, ]
    reg <- lm( val ~ p, data = p___ )
    slope <- c( slope, round( coef(reg)[2], 2 ) )
  }
  Slope$Kapur_norm <- cbind( RQ_, slope )
  Slope$Kapur_norm <- reshape( Slope$Kapur_norm, idvar = "AR", timevar = "BQ", direction = "wide" )
  
  p_$Q <- as.factor( p_$Q )
  
  fig <- plot_ly( x = p_$p, y = p_$R, z = p_$val, color = p_$Q )
  fig <- add_markers( fig )
  fig <- layout( fig, title = list( text = paste0( Dist, " normalized distribution" ), y = 0.8 ), 
                 scene = list(xaxis = list(title = 'Probability', range=c( 0, 1 ) ),
                              yaxis = list(title = 'A'),
                              zaxis = list(title = 'Kapur' ),
                              camera = list( eye = list( x = 1.8, y = -1.5, z = 0.5 ) )
                 ),
                 legend = list( x = 1, y = 0.1 ), 
                 margin = list( l = 0, r = 0, b = 0, t = 0, pad = 0 ),
                 annotations = list( x = 1, y = 0.25, text = "B", xref = 'paper', yref = 'paper', showarrow = FALSE )
  )
  
  FIG[[ 16 ]] <- fig
  
  
  #### Out ####
  names( FIG ) <- c( "Hist_dist","Hist_prob","Hist_dist_norm","Hist_prob_norm",
                     "Shannon_Miss_Gini","Renyi","Tsallis","Sharma-Mittal3D","Sharma-Taneja3D","Kapur3D",
                     "Shannon_Miss_Gini_norm","Renyi_norm","Tsallis_norm", 
                     "Sharma-Mittal3D_norm","Sharma-Taneja3D_norm","Kapur3D_norm"
                    )
  
  return( c( FIG, Slope ) )
  
}

CrossValid <- function( Y_name, X_names, dat, type, depth, min_obs, Qval, Alpha, Beta, weights, cost, overfit, cp, cf, k_fold, seed ){
  
  Shannon <- RenyiTsallis <- MittalTanejaKapur <- NULL
  
  if( "Shannon" %in% type ){
    
    if( c("prune") %in% overfit ){
      
      Shannon1 <- expand.grid( Type = "Shannon", Q = 1, Alpha = 1, Beta = 1, Depth = depth, Min_obs = min_obs, 
                               Overfit = overfit[which("prune"!=overfit)], Cp = cp, Cf = 0.25, stringsAsFactors = F )
      Shannon2 <- expand.grid( Type = "Shannon", Q = 1, Alpha = 1, Beta = 1, Depth = depth, Min_obs = min_obs, 
                               Overfit = "prune", Cp = 0, Cf = cf, stringsAsFactors = F )
      Shannon <- rbind( Shannon1, Shannon2 )
      
    }else{
      
      Shannon <- expand.grid( Type = "Shannon", Q = 1, Alpha = 1, Beta = 1, Depth = depth, Min_obs = min_obs, 
                              Overfit = overfit, Cp = cp, Cf = 0.25, stringsAsFactors = F )
      
    }
    
  }
  if( any(c("Renyi","Tsallis") %in% type) ){
    
    if( c("prune") %in% overfit ){
      
      RenyiTsallis1 <- expand.grid( Type = type[type%in%c("Renyi","Tsallis")], Q = Qval, Alpha = 1, Beta = 1, 
                                    Depth = depth, Min_obs = min_obs, Overfit = overfit[which("prune"!=overfit)], 
                                    Cp = cp, Cf = 0.25, stringsAsFactors = F )
      RenyiTsallis2 <- expand.grid( Type = type[type%in%c("Renyi","Tsallis")], Q = Qval, Alpha = 1, Beta = 1, 
                                    Depth = depth, Min_obs = min_obs, Overfit = "prune", Cp = 0, Cf = cf, 
                                    stringsAsFactors = F )
      RenyiTsallis <- rbind( RenyiTsallis1, RenyiTsallis2 )
      
    }else{
      
      RenyiTsallis <- expand.grid( Type = type[type%in%c("Renyi","Tsallis")], Q = Qval, Alpha = 1, Beta = 1,
                                   Depth = depth, Min_obs = min_obs, Overfit = overfit, Cp = cp, Cf = 0.25, stringsAsFactors = F )
      
    }
    
  }
  if( any(c("Sharma-Mittal", "Sharma-Taneja", "Kapur") %in% type) ){
    
    if( c("prune") %in% overfit ){
      
      MittalTanejaKapur1 <- expand.grid( Type = type[type%in%c("Sharma-Mittal","Sharma-Taneja","Kapur")],
                                         Q = 1, Alpha = Alpha, Beta = Beta, Depth = depth, Min_obs = min_obs, 
                                         Overfit = overfit[which("prune"!=overfit)], Cp = cp, Cf = 0.25, stringsAsFactors = F )
      MittalTanejaKapur2 <- expand.grid( Type = type[type%in%c("Sharma-Mittal","Sharma-Taneja","Kapur")], 
                                         Q = 1, Alpha = Alpha, Beta = Beta, Depth = depth, Min_obs = min_obs, 
                                         Overfit = "prune", Cp = 0, Cf = cf, stringsAsFactors = F )
      MittalTanejaKapur <- rbind( MittalTanejaKapur1, MittalTanejaKapur2 )
      
    }else{
      
      MittalTanejaKapur <- expand.grid( Type = type[type%in%c("Sharma-Mittal","Sharma-Taneja","Kapur")], 
                                        Q = 1, Alpha = Alpha, Beta = Beta, Depth = depth, Min_obs = min_obs, 
                                        Overfit = overfit, Cp = cp, Cf = 0.25, stringsAsFactors = F )
      
    }
    
  }
  exclude <- (MittalTanejaKapur$Type == "Kapur" & ( (MittalTanejaKapur$Alpha == 1) | (MittalTanejaKapur$Alpha <= 0) | 
                                                     (MittalTanejaKapur$Beta <= 0) | (MittalTanejaKapur$Alpha + MittalTanejaKapur$Beta - 1 <= 0) ) ) | 
    (MittalTanejaKapur$Type == "Sharma-Taneja" & MittalTanejaKapur$Alpha == MittalTanejaKapur$Beta) |
    (MittalTanejaKapur$Type == "Sharma-Mittal" & ( (MittalTanejaKapur$Alpha == MittalTanejaKapur$Beta) | 
                                                     MittalTanejaKapur$Beta == 1 ) )
  MittalTanejaKapur <- MittalTanejaKapur[!exclude,]
  
  Grid <- rbind( Shannon, RenyiTsallis, MittalTanejaKapur )
  
  set.seed( seed )
  
  cross_validation <- createFolds( dat[,Y_name], k_fold, T )
  cross_validation <- lapply( cross_validation, function( valid, tar ){ 
    temp <- rep(1, length(tar)); temp[valid] <- 2; temp 
  }, tar = dat[,Y_name])
  
  train <- data.frame( matrix( 0, length(cross_validation) * nrow(Grid), ncol(Grid) + 5 ) )
  colnames(train) <- c( "Cross", colnames(Grid), "Accuracy", "Kappa", "Nleaves", "Nclass" )
  valid <- train
  
  k <- 1
  for( i in 1:length(cross_validation) ){
    
    for( j in 1:nrow(Grid) ){
      
      Tree <- ImbTreeEntropy( Y_name = Y_name, X_names = X_names, data = dat[ cross_validation[[i]] == 1, ], type = Grid[j,"Type"], 
                              depth = Grid[j,"Depth"], min_obs = Grid[j,"Min_obs"], cp = Grid[j,"Cp"], cf = Grid[j,"Cf"], 
                              overfit = Grid[j,"Overfit"], weights = weights, cost = cost, class_th = "equal",
                              entropy_par = if(Grid[j,"Type"] %in% c("Sharma-Mittal", "Sharma-Taneja", "Kapur")){c(Grid[j,"Alpha"], Grid[j,"Beta"])}else{Grid[j,"Q"]} )
      
      pred_train <- PredictTree( Tree, dat[ cross_validation[[i]] == 1, ] )
      pred_valid <- PredictTree( Tree, dat[ cross_validation[[i]] == 2, ] )
      
      conf_train <- confusionMatrix( pred_train$Class, dat[ cross_validation[[i]] == 1, Y_name ] )$overall[1:2]
      conf_valid <- confusionMatrix( pred_valid$Class, dat[ cross_validation[[i]] == 2, Y_name ] )$overall[1:2]
      
      train[k, "Cross"] <- i
      valid[k, "Cross"] <- i
      
      train[k, !colnames(train) %in% c("Cross","Accuracy","Kappa") ] <- Grid[j,]
      valid[k, !colnames(valid) %in% c("Cross","Accuracy","Kappa")] <- Grid[j,]
      
      train[k, c("Accuracy", "Kappa")] <- conf_train[c("Accuracy", "Kappa")]
      valid[k, c("Accuracy", "Kappa")] <- conf_valid[c("Accuracy", "Kappa")]
      
      valid[k, "Nleaves"] <- train[k, "Nleaves"] <- length(Tree$Get("pathString", filterFun = isLeaf))
      valid[k, "Nclass"] <- train[k, "Nclass"] <- length(unique(Tree$Get("Class", filterFun = isLeaf)))
      
      print( sprintf("%s out of %s", k, nrow(train) ) )
      k <- k + 1
      
    }
    
  }
  
  train[, "Nleaves"] <- as.numeric(train[, "Nleaves"])
  valid[, "Nleaves"] <- as.numeric(valid[, "Nleaves"])
  
  train_agg <- aggregate( train[,c("Accuracy","Kappa","Nleaves","Nclass")], as.list(train[,2:10]), mean )
  valid_agg <- aggregate( valid[,c("Accuracy","Kappa","Nleaves","Nclass")], as.list(valid[,2:10]), mean )
  
  return( list(Train = train, TrainAgg = train_agg, Valid = valid, ValidAgg = valid_agg) )
  
}

sensPlot <- function( dat, type = "Shannon", cond = NULL, title = "" ){
  
  cond1 <- paste( sprintf( 'dat[,"Type"] == "%s" &', type ), cond )
  condS <- paste( 'dat[,"Type"] == "Shannon" &', cond )
  
  datS <- dat[ eval( parse( text = condS ) ), ]
  dat <- dat[ eval( parse( text = cond1 ) ), ]

  if( type %in% c( "Renyi", "Tsallis" ) ){
    
    fig <- plot_ly( x = dat[,"Q"], y = dat[,"Accuracy"], type = 'scatter', 
                    mode = 'lines+markers', name = paste0( "Accuracy: ", type ) )
    fig <- add_trace( fig, x = dat[,"Q"], y = dat[,"Nleaves"], name = paste0( "Leaves: ", type ), 
                      yaxis = "y2", mode = "lines", type = "scatter", line = list( dash = 'dash' ) )
    
    fig <- add_trace( fig, x = dat[,"Q"], y = datS[,"Accuracy"], name = "Accuracy: Shannon", 
                      yaxis = "y", mode = "lines+markers", type = "scatter" )
    fig <- add_trace( fig, x = dat[,"Q"], y = datS[,"Nleaves"], name = "Leaves: Shannon", 
                      yaxis = "y2", mode = "lines", type = "scatter", line = list( dash = 'dash' ) )
    
  }else{
    
    uniq <- unique( dat[,"Beta"] )
    
    if( type == "Sharma-Mittal" ){
      name <- c( "R", "Q" )  
      uniq <- uniq[ seq( 1, length(uniq), by = 2 ) ]
    }else if( type == "Sharma-Taneja" ){
      name <- c( "A", "B" )
      uniq <- uniq[-1]
      uniq <- uniq[ seq( 1, length(uniq), by = 2 ) ]
    }else{
      name <- c( "A", "B" )
      uniq <- uniq[ seq( 1, length(uniq), by = 2 ) ]
    }
    
    for( i in 1:length(uniq) ){
      
      datS <- dat[ dat[,"Beta"] == uniq[i], ]
      if( i == 1 ){
        
        fig <- plot_ly( x = datS[,"Alpha"], y = datS[,"Accuracy"], type = 'scatter', 
                        mode = 'lines+markers', name = paste0( "Accuracy: ", name[2], " ", uniq[i] ) )
        
      }else{
        
        fig <- add_trace( fig, x = datS[,"Alpha"], y = datS[,"Accuracy"], name = paste0( "Accuracy: ", name[2], " ", uniq[i] ), 
                          yaxis = "y", mode = "lines+markers", type = "scatter" )
        
      }
      
      fig <- add_trace( fig, x = datS[,"Alpha"], y = datS[,"Nleaves"], name = paste0( "Leaves: ", name[2], " ", uniq[i] ), 
                        yaxis = "y2", mode = "lines", type = "scatter", line = list( dash = 'dash' ) )
      
    }
    
  }

  
  fig <- layout( fig, title = title, 
                 xaxis = list( title = name[1], range=c( 0, 5 ) ), 
                 yaxis = list( title = "Accuracy" ),
                 yaxis2 = list( overlaying = "y", side = "right", title = "Number of leaves" ) )

  return( fig )
  
}

##### Plots one variable ####

## Normal ( 0, 1 )
set.seed(123)

distTypeNormal  <- list( Normal = list( "rnorm", NULL ) )
outNormal <- SimVar( n = 1000, distType = distTypeNormal )

ppNormal <- Plots( outNormal, T, T )

## Beta ( alpha, beta )
set.seed(123)

distTypeBeta  <- list( Beta = list( "rbeta", c( 0.5,0.5 ) ) )
distTypeBeta2  <- list( Beta = list( "rbeta", c( 5, 1 ) ) )

outBeta <- SimVar( n = 1000, distType = distTypeBeta )
outBeta2 <- SimVar( n = 1000, distType = distTypeBeta2 )

ppBeta  <- Plots( outBeta, T, T )
ppBeta2  <- Plots( outBeta2, T, T )

## Cauchy ( location = 0, scale = 1 )
set.seed(123)

distTypeCauchy <- list( Cauchy = list( "rcauchy", NULL ) )
outCauchy <- SimVar( n = 1000, distType = distTypeCauchy )

ppCauchy <- Plots( outCauchy, T, T )

## Exponential ( lambda = 1 )
set.seed(123)

distTypeExponential <- list( Exponential = list( "rexp", NULL ) )
outExp <- SimVar( n = 1000, distType = distTypeExponential )

ppExp <- Plots( outExp, T, T )

## Uniform
set.seed(123)

distTypeUniform <- list( Uniform = list( "runif", NULL ) )
outUniform <- SimVar( n = 1000, distType = distTypeUniform )

ppUniform <- Plots( outUniform, T, T )

#### Interactive learning ####
set.seed(123)
distLinearCombination <- list( Normal = list( "rnorm", NULL ), Beta = list( "rbeta", c(0.5,0.5) ),
                               Beta2 = list( "rbeta", c( 5, 1 ) ), Cauchy = list( "rcauchy", NULL ), 
                               Exponential = list( "rexp", NULL ), Uniform = list( "runif", NULL ) )
linearCombination <- SimVar( n = 1000, distType = distLinearCombination )

Data <- data.frame( linearCombination$Data, 
                    Y = factor( ifelse( linearCombination$ProbsNorm <= 0.5, "_0_", "_1_" ) ) )

write.table( Data, "Data.txt", sep = ";", row.names = F )

library( "remotes" )
install_github( "KrzyGajow/ImbTreeEntropy" )

library( "ImbTreeEntropy" )
runShinyImbTreeEntropy()

#### Sensitivity Analysis ####

SensAna <-  CrossValid( Y_name = "Y", X_names = colnames(Data)[1:6], dat = Data, 
                        type = c( "Shannon", "Renyi", "Tsallis", "Sharma-Mittal", "Sharma-Taneja", "Kapur" ), 
                        depth = seq( 5, 10 ), 
                        min_obs = c( 10, 50, 100 ), 
                        Qval = seq( 0, 5, 0.5 ), 
                        Alpha = seq( 0, 5, 0.5 ), 
                        Beta = seq( 0, 5, 0.5 ),
                        weights = NULL, 
                        cost = NULL,
                        overfit = "leafcut", 
                        cp = seq( 0, 0, 0 ), 
                        cf = seq( 0, 0, 0 ),
                        k_fold = 5, 
                        seed = 666 )

sensPlot( dat = SensAna$ValidAgg, 
          type = "Renyi",
          cond = 'dat[,"Depth"] == 5 & dat[,"Min_obs"] == 50',
          title = "Renyi" )

sensPlot( dat = SensAna$ValidAgg, 
          type = "Tsallis",
          cond = 'dat[,"Depth"] == 5 & dat[,"Min_obs"] == 50',
          title = "Tsallis" )

sensPlot( dat = SensAna2$ValidAgg, 
          type = "Sharma-Mittal",
          cond = 'dat[,"Depth"] == 5 & dat[,"Min_obs"] == 50',
          title = "Sharma-Mittal" )

sensPlot( dat = SensAna2$ValidAgg, 
          type = "Sharma-Taneja",
          cond = 'dat[,"Depth"] == 5 & dat[,"Min_obs"] == 50',
          title = "Sharma-Taneja" )

sensPlot( dat = SensAna2$ValidAgg, 
          type = "Kapur",
          cond = 'dat[,"Depth"] == 5 & dat[,"Min_obs"] == 50',
          title = "Kapur" )
