predict0<- function( mod, meig0, x0 = NULL ){
    if( mod$other$model == "esf" ){
    	meig0$sf	<- as.matrix( meig0$sf[, mod$other$sf_id ] )
    }

    if( is.null( dim( mod$r ) ) ){
    	sf_pred	<- meig0$sf * mod$r[  1 ]
    } else {
    	sf_pred	<- meig0$sf %*% mod$r[, 1 ]
    }

    if( is.null( mod$other$x_id ) ){
    	xb_pred	<- c( as.matrix( mod$b[ 1 ] ) )
    	pred	<- xb_pred + sf_pred
    } else {
    	if( is.null( x0 ) ){
    		message( " Note: Only spatial component (sf) is interpolated because x0 is missing")
    		pred	<- sf_pred
   	 } else {
   		x02	<- as.matrix( x0 )
  		if( is.numeric( x02 ) == F ){
    			mode( x02 ) <- "numeric"
    		}

   	 	xb_pred	<- as.matrix( cbind( 1, x02[ ,mod$other$x_id ] ) ) %*% mod$b[, 1 ]
    		pred	<- xb_pred + sf_pred
    		pred	<- as.data.frame( cbind( pred, xb_pred, sf_pred ) )
    		names( pred )<- c( "pred", "xb", "sf" )
    	}
    }
    return( pred )
}

