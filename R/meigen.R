meigen		<- function( coords, threshold = 0, enum = NULL, cmat = NULL, longlat = FALSE){
    if( threshold > 1 | threshold < 0 ){
    	stop( "threshold must lie between 0 and 1" )
    } else if ( threshold ==1) {
    	threshold <- threshold - 1e-07
    }
    if( is.null( cmat )){
    	if( longlat == TRUE){
    		D	<- rdist.earth( coords )
    	} else {
    		D	<- rdist( coords )
    	}
    	h	<- max( spantree( D )$dist )
    	C	<- exp( -D / h )
    } else {
    	if( isSymmetric( cmat )==F ) {
    		C 	<- 0.5 * ( cmat + t( cmat ) )
    		message( " Note:" )
    		message( "   cmat is symmetrized by ( cmat + t( cmat ) ) / 2" )
    	}
    	C	  <- as.matrix( cmat )
    	diag( C ) <- 0
	cmat	  <- NULL
    }
    diag( C )	<- 0
    Cmean	<- apply( C, 1, mean )
    MCM		<- t( C - Cmean ) - Cmean + mean( Cmean )
    eigenC	<- eigen( MCM )
    eigenC$values	<- Re( eigenC$values )
    eigenC$vectors	<- Re( eigenC$vectors )
    sel		<- ( eigenC$values / max( eigenC$values ) >= threshold + 1e-07 )
    if( is.null( enum ) == F ){
    	if( sum(sel) > enum ){
    		sel[ -c(1:enum) ] <- FALSE
    	}
    }
    sf		<- as.matrix( eigenC$vectors[ , sel ] )
    ev		<- eigenC$values [ sel ]

    mes		<- paste( " ", length( ev ), "/", length( sf[ , 1 ] ), " eigen-pairs are extracted", sep = "")
    message( mes )
    return( list( sf = sf, ev = ev, ev_full = eigenC$values, fast = 0 ) )
}

