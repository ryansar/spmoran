resf  	<- function( y, x = NULL, meig, method = "reml" ){

    lik_resf	<- function( par0, ev, M, m, yy, n, nx, ne, emet ){
    	par	<- par0 ^ 2
    	evv	<- ev ^ par[ 1 ] * ( sum( ev ) / sum( ev ^ par[ 1 ] ) )
    	evSqrt	<- par[ 2 ] * sqrt( evv )
    	Mw	<- t( M[ -( 1:nx ), -( 1:nx ) ] * evSqrt ) * evSqrt
    	M[  -( 1:nx ), -( 1:nx ) ]	<- Mw + diag( ne )
    	M[  -( 1:nx ),    1:nx   ]	<- M[ -( 1:nx ), 1:nx ] * evSqrt
    	M[     1:nx  , -( 1:nx ) ]	<- t( M[ -( 1:nx ), 1:nx ] )
    	M0				<- M
    	M0[ -( 1:nx ), -( 1:nx ) ]	<- Mw
    	m[  -( 1:nx) ]			<- m[ -( 1:nx ) ] * evSqrt

    	test    <-try( Minv	<- solve( M, tol = 1e-25 ) )
    	if( class( test ) == "try-error" ){
    		loglik  <- Inf
    	} else {
    		b	<- Minv %*% m
    		sse	<- yy - 2 * t( b ) %*% m + t( b ) %*% M0 %*% b
    		dd	<- sse + sum( b[ -( 1:nx ) ] ^ 2 )
    		if( emet == "reml" ){
    			term1	<- determinant( M )$modulus
    			term2	<- ( n - nx ) * ( 1 + log( 2 * pi * dd / ( n - nx ) ) )
    		} else if( emet == "ml" ){
    			term1	<- determinant( as.matrix( M[ -( 1:nx ), -( 1:nx ) ] ) )$modulus
    			term2	<- n * ( 1 + log( 2 * pi * dd / n ) )
    		}
    		loglik	<- term1 + term2
    	}
    	return( loglik )
    }

    n		<- length( y )
    if( is.null( x ) ){
    	X	<- as.matrix( rep( 1, n ) )
    	xname	<- "(Intercept)"
    	x_id	<- NULL
    } else {
    	X00	<- as.matrix( x )
    	if( is.numeric( X00 ) == F ){
    		mode( X00 ) <- "numeric"
    	}
    	x_id	<- apply( X00, 2, sd ) != 0
    	if( sum( x_id ) == 0 ){
    		X	<- as.matrix( rep( 1, n ) )
    		xname	<- "(Intercept)"
    		x_id	<- NULL
    	} else {
    		X0	<- X00[ , x_id ]
    		X	<- as.matrix( cbind( 1, X0 ) )
    		xname	<- c( "(Intercept)", names( as.data.frame( X0 ) ) )
    	}
    }

    ev		<- meig$ev
    nx		<- dim( X )[ 2 ]
    ne		<- length( ev )
    yy     	<- sum( y ^ 2 )
    XX		<- crossprod(  X )
    Xy		<- crossprod(  X, y )
    EX		<- crossprod( meig$sf, X )
    Ey		<- crossprod( meig$sf, y )
    if( meig$other$fast == 0 ){
    	EE	<- diag( ne )
    } else if( meig$other$fast == 1){
    	EE	<- crossprod( meig$sf )
    }
    M		<- as.matrix( rbind( cbind( XX, t( EX ) ), cbind( EX, EE ) ) )
    m		<- c( Xy, Ey )
    res		<- optim( fn = lik_resf, c( 1, 1 ), ev = ev, M = M, m = m, yy = yy,
    			   n = n, nx = nx, ne = ne, emet = method )
    par		<- res$par ^ 2
    loglik	<- ( -1 / 2 ) * res$value
    evv		<- ev ^ par[ 1 ] * ( sum( ev ) / sum( ev ^ par[ 1 ] ) )
    evSqrt	<- par[ 2 ] * sqrt( evv )
    sf2		<- t( t( meig$sf ) * evSqrt )
    X2		<- as.matrix( cbind( X, sf2 ) )
    XX2		<- crossprod( X2 )
    diag( XX2 )[ -( 1:nx ) ] <- diag( XX2 )[ -( 1:nx ) ] + 1
    XXinv2	<- solve( XX2 )
    XXinv2_X	<- XXinv2 %*% t( X2 )
    b		<- XXinv2_X %*% y
    b[ -( 1:nx ) ] <- b[ -( 1:nx ) ] * evSqrt
    pred	<- as.matrix( cbind( X, meig$sf ) ) %*% b
    resid	<- y - pred
    SSE		<- sum( resid ^ 2 )
    SSY		<- sum( ( y - mean( y ) ) ^ 2 )
    sig		<- SSE / ( n - nx )
    bse		<- sqrt( sig ) * sqrt( diag( XXinv2 ) )
    SF		<- meig$sf %*% b[ -( 1:nx ) ]

    np		<- nx + 1 + 3
    AIC		<- -2 * loglik + np * 2
    BIC		<- -2 * loglik + np * log( n )
    r2_0	<- 1 - SSE / SSY
    r2		<- 1 - ( 1- r2_0 ) * ( n - 1 ) / ( n - np - 1)

    bt		<- b[ 1:nx ] / bse[ 1:nx ]
    df		<- sum(t(X2)*XXinv2_X)
    bp		<- 2 - 2 * pt( abs( bt ), df = n - df )

    b_par	<- data.frame( Estimate = b[ 1:nx ], SE = bse[ 1:nx ], t_value = bt, p_value = bp )
    rownames( b_par ) <- xname

    r_par		<- data.frame( b[ -( 1:nx ) ] )
    names( r_par )	<- "Estimate"
    rownames( r_par )	<- paste( "r", 1:ne, sep = "" )

    par[ 2 ]	<- par[ 2 ] * sqrt( sig )
    sf_par	<- data.frame( par = par[ c( 2, 1 ) ] )
    rownames( sf_par ) <- c( "shrink_sf_SE", "shrink_sf_alpha" )

    e_stat	<- data.frame( stat = c( sqrt( sig ), r2, loglik, AIC, BIC ) )
    if( method == "reml" ){
    	rownames( e_stat ) <- c( "resid_SE", "adjR2(cond)", "rlogLik", "AIC", "BIC" )
    } else if( method == "ml" ){
    	rownames( e_stat ) <- c( "resid_SE", "adjR2(cond)", "logLik", "AIC", "BIC" )
    }
    other	<- list( x_id = x_id, model = "resf", par0 = res$par )

    return( list( b = b_par, s = sf_par, e = e_stat,
    		  r = r_par, sf = SF, pred = pred, resid = resid, other = other ) )
}
