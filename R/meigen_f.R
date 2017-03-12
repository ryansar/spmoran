meigen_f	<- function( coords, enum = 200 ){
    n		<- dim( coords )[ 1 ]
    if( enum < 1 | enum != floor( enum ) ){
    	stop( " enum must be a positive integer smaller than n" )
    }
    if( enum >= n ){
    	stop( " enum must be less than n" )
    }


    if( enum >= n - 1 ){
    	meig	<- meigen(coords=coords)
    	ev_ap	<- meig$ev
    	ev_ap0	<- meig$ev_full
    	sf_ap	<- meig$sf
    	fast	<- 0
    } else {
    	coordk	<- kmeans( coords, centers = enum + 1 )$centers
    	D	<- rdist( coordk )
    	h	<- max( spantree( D )$dist )
    	C	<- exp( -D / h )
    	Cmean	<- apply( C, 1, mean )
    	MCM	<- t( C - Cmean ) - Cmean + mean( Cmean )
    	eigenC	<- eigen( MCM )
    	sf0	<- eigenC$vectors[ , eigenC$values > 1e-08 ]
    	ev0	<- eigenC$values [   eigenC$values > 1e-08 ]
    	ev_ap0	<- ev0 * ( n / enum ) - 1
    	sf_ap0	<- sf0 %*% diag( 1 / ev0 )
    	sf_ap0	<- t( exp( -rdist( coordk, coords ) / h ) - Cmean ) %*% sf_ap0
    	sf_ap0	<- t( t( sf_ap0 ) - colMeans( sf_ap0 ) )
    	ev_ap	<- ev_ap0[   ev_ap0 > 1e-08 ]
    	sf_ap	<- as.matrix( sf_ap0[ , ev_ap0 > 1e-08 ] )
    	mes	<- paste( " ", length( ev_ap ), "/", length( sf_ap[ , 1 ] ), " eigen-pairs are approximated", sep = "" )
    	message( mes )
    	fast	<- 1
    }
    return( list( sf  =sf_ap, ev = ev_ap, ev_full = ev_ap0, fast = fast ) )
}



