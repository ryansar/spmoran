
meigen0	<- function( meig, coords0 ){

    if( meig$other$model == "other" ){
    	stop( "meigen0 is not supported for user-specified spatial proximity matrix" )
    }

    if( meig$other$fast == 0){
    	sfk	<- meig$sf
    	evk	<- meig$ev
    	coordk	<- meig$other$coords
    } else {
    	sfk	<- meig$other$sfk
    	evk	<- meig$other$evk
    	coordk	<- meig$other$coordk
    }
    h		<- meig$other$h
    nm		<- dim( coords0 )[ 1 ]
    sfk		<- sfk %*% diag( 1 / ( evk + 1 ) )
    Dk		<- rdist( coordk, coords0 )
    	if( meig$other$model == "exp" ){
    		sfk	<- t( exp( -Dk / h ) - meig$other$Cmean ) %*% sfk
    	} else if( meig$other$model == "gau" ){
    		sfk	<- t( exp( - ( Dk / h ) ^ 2 ) - meig$other$Cmean ) %*% sfk
    	} else if( meig$other$model == "sph" ){
    		sfk	<- t( ifelse( Dk < h , 1 - 1.5 * ( Dk / h ) + 0.5 * ( Dk / h ) ^ 3, 0 ) - meig$other$Cmean ) %*% sfk
    	}
return( list( sf = sfk, ev = meig$ev, ev_full = meig$ev_full )  )
}
