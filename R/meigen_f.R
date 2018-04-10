meigen_f	<- function( coords, model = "exp", enum = 200 ){			
    n		<- dim( coords )[ 1 ]		
    if( enum < 1 | enum != floor( enum ) ){				
    	stop( " enum must be a positive integer" )			
    }				
    				
    if( enum >= n - 1 ){				
    	result	<- meigen(coords = coords, model = model )
    	result$other$coordk	<- NULL
    	result$other$sfk	<- NULL
    	result$other$evk	<- NULL

    } else {				
    	coordk	<- kmeans( coords, centers = enum + 1 )$centers		
    	D	<- rdist( coordk )		
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
    				
    	Cmean	<- apply( C, 1, mean )		
    	MCM	<- t( C - Cmean ) - Cmean + mean( Cmean )		
    	eigenC	<- eigen( MCM )		
    	sf0	<- eigenC$vectors[ , eigenC$values > 1e-08 ]		
    	ev0	<- eigenC$values [   eigenC$values > 1e-08 ]		
    	ev_full	<- ev0 * ( n / enum ) - 1		
    	ev_ap	<- ev_full[   ev_full > 1e-08 ]		

    	sf_ap0	<- sf0 %*% diag( 1 / ev0 )		
    	sf_ap0	<- t( exp( -rdist( coordk, coords ) / h ) - Cmean ) %*% sf_ap0		
    	sf_ap0	<- t( t( sf_ap0 ) - colMeans( sf_ap0 ) )		
    	sf_ap	<- as.matrix( sf_ap0[ , ev_full > 1e-08 ] )		
   	sf0	<- as.matrix( sf0   [ , ev_full > 1e-08 ] )
    	ev0	<- ev0   [   ev_full > 1e-08 ] -1
    	other	<- list( coordk = coordk, sfk = sf0, evk = ev0, Cmean = Cmean, h = h, model = model, fast = 1 )
    	mes	<- paste( " ", length( ev_ap ), "/", length( sf_ap[ , 1 ] ), " eigen-pairs are approximated", sep = "" )		
    	message( mes )			
    	result	<- list( sf  =sf_ap, ev = ev_ap, ev_full = ev_full, other = other )
    }				
    return( result )				
}				
