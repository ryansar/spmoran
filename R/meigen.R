meigen		<- function( coords, model = "exp", threshold = 0, enum = NULL, cmat = NULL ){		
    if( threshold > 1 | threshold < 0 ){				
    	stop( "threshold must lie between 0 and 1" )			
    } else if ( threshold ==1) {				
    	threshold <- threshold - 1e-07			
    }				
    if( is.null( cmat )){				
    	D	<- rdist( coords )		
    	h	<- max( spantree( D )$dist )		
    	if( model == "exp" ){			
    		C	<- exp( -D / h )	
    	} else if( model == "gau" ){			
    		C	<- exp( -( D / h) ^ 2 )	
    	} else if( model == "sph" ){			
    		C	<- ifelse( D < h , 1 - 1.5 * (D / h ) + 0.5 * ( D / h ) ^ 3, 0 )	
    	} else {
    		stop( "model is not specified appropriately" )
    	}		
    } else {				
    	if( isSymmetric( cmat )==F ) {			
    		C 	<- 0.5 * ( cmat + t( cmat ) )	
    		message( " Note:" )		
    		message( "   cmat is symmetrized by ( cmat + t( cmat ) ) / 2" )		
    	}			
    	C	<- as.matrix( cmat )	
    	model	<- "other"
    	coords	<- NULL
    	h	<- NULL
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
    other	<- list( coords = coords, Cmean = Cmean, h = h, model = model, fast = 0 )
    mes		<- paste( " ", length( ev ), "/", length( sf[ , 1 ] ), " eigen-pairs are extracted", sep = "" )		
    message( mes )
    return( list( sf = sf, ev = ev, ev_full = eigenC$values, other = other ) )				
}

