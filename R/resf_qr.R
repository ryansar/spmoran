resf_qr <- function( y, x = NULL, meig, tau = NULL, boot = TRUE, iter = 200 ){

    resf_boot  	<-function( tau_sel, q_sel, y, X, M, meig, mod, iter ){

      re_esfb	<-function( y, X, XSF0, M, meig ){

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
            if( ( emet == "reml" ) | ( emet == "pls" ) ){
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

        ev	<- meig$ev
        n	<- length( y )
        nx	<- dim( X )[ 2 ]
        ne	<- length( ev )
        yy	<- sum( y^2 )
        m	<- c( crossprod( X, y ), crossprod( meig$sf, y ) )
        res 	<- optim( fn = lik_resf, c( 1, 1 ), ev = ev, M = M, m = m,
                       yy = yy, n = n, nx = nx, ne = ne, emet = "pls" )
        par	<- res$par ^ 2
        loglik 	<- ( -1 / 2 ) * res$value
        evv     <- ev ^ par[ 1 ] * ( sum( ev ) / sum( ev ^ par[ 1 ] ) )
        evSqrt  <- par[ 2 ] * sqrt( evv )
        sf2     <- t( t( meig$sf ) * evSqrt )
        X2     	<- as.matrix( cbind( X, sf2 ) )
        XX2  	<- crossprod( X2 )
        diag( XX2 )[ -( 1:nx ) ] <- diag( XX2 )[ -( 1:nx ) ] + 1
        XXinv2	<- solve( XX2 )
        b	<- XXinv2 %*% t( X2 ) %*% y
        b[ -( 1:nx ) ]	<- b[ -( 1:nx ) ] * evSqrt
        sig		<- sum( ( y - XSF0 %*% b ) ^ 2 ) / ( n - nx )
        par[ 2 ]	<- par[ 2 ] * sqrt( sig )
        return(list( b = b[ 1:nx ], par = par ) )
      }

      b	<- mod$b[  , 1 ]
    	s	<- mod$s[ 2:1, ]
    	sd	<- mod$e[ 1, 1 ]
    	ev	<- meig$ev
    	n	<- length( y )
    	ne	<- length( ev )
    	nx	<- dim( X )[ 2 ]

    	evv	<- ev ^ s[ 1 ] * ( sum( ev ) / sum( ev ^ s[ 1 ] ) )
    	evSqrt	<- s[ 2 ] * sqrt( evv )
    	sf2	<- t( t( meig$sf ) * evSqrt )
    	XSF0	<- as.matrix( cbind( X, meig$sf ) )
    	f0	<- density( y )
    	fq0	<- approx( f0$x, f0$y, q_sel )$y

    	res<- foreach( i = 1:iter, .combine = "rbind" ) %dopar% {
    		usim	<- rnorm( n , sd = sd )
    		rsim	<- rnorm( ne, sd = 1  )
    		y_b	<- X %*% b + sf2 %*% rsim + usim
   		f_b	<- density( sample(y, replace=TRUE) )
    		fq_b	<- approx( f_b$x, f_b$y, q_sel )$y
    		RIF_b	<- fq0 / fq_b * ( y_b - q_sel) + q_sel
    		sfres_sim	<- re_esfb( y = RIF_b, X = X, XSF0 = XSF0, M = M, meig = meig )
    		c( sfres_sim$b, sfres_sim$par[ 2:1 ] )
    	}
    	B	<- as.matrix( res[ , 1:nx ] )
    	S	<- as.matrix( res[ ,( nx + 1): (nx + 2) ] )
    	return( list( B = B, S = S ) )
    }

    indic <- function( y, q_sel, tau_sel ){
    	if( tau_sel == 1 ){
    		res <-ifelse( y < q_sel, 1, 0 )
    	} else {
    		res <-ifelse( y <= q_sel, 1, 0 )
    	}
    	return(res)
    }

    if(is.null( tau ) ) tau <- 1:9 / 10
    q	  <- quantile( y, probs = tau )
    f	  <- density( y )
    fq	<- approx( f$x, f$y, q )$y
    n   <- length( y )
    ne  <- length( meig$ev )

    if( is.null( x ) ){
    	X 	<- as.matrix( rep( 1, n ) )
    	xname 	<- "(Intercept)"
    } else {
    	X00	<- as.matrix( x )
    	if( is.numeric( X00 ) == F ){
    		mode( X00 ) <- "numeric"
    	}
    	xind	<- apply( X00, 2, sd ) != 0
    	X0	<- X00[ , xind ]
    	X	<- as.matrix( cbind( 1, X0 ) )
        if( sum( xind ) == 0 ){
    		xname <- "(Intercept)"
        } else {
    		xname <- c( "(Intercept)", names( as.data.frame( X0 ) ) )
    	}
    }

    if( boot == TRUE ){
    	XX	<- crossprod( X )
    	EX	<- crossprod( meig$sf, X )
    	if( meig$other$fast == 0 ){
    		EE	<- diag( length( meig$ev ) )
    	} else if( meig$other$fast == 1 ){
    		EE	<- crossprod( meig$sf )
    	}
    	M	<- as.matrix( rbind( cbind( XX, t( EX ) ), cbind( EX, EE ) ) )
    	cl	<- makeCluster( detectCores() )
    	registerDoParallel( cl )
    }

    SFb  	<- NULL
    SFr  	<- NULL
    SFs		<- NULL
    SFe		<- NULL
    SFb_boot	<-list( NULL )
    SFs_boot	<-list( NULL )
    SFb_boot2 <-list( NULL )
    SFs_boot2 <-list( NULL )
    probs	    <- c( 0.025, 0.975 )
    for( j in 1:length( tau ) ){
      RIF	<- q[ j ] + ( tau[ j ] - indic( y = y, q_sel = q[ j ], tau_sel = tau[ j ] ) ) / fq[ j ]
    	mod	<- resf( y = RIF, meig = meig, x = X )
    	SFb	<- cbind( SFb	, mod$b[ , 1 ] )
    	SFr	<- cbind( SFr	, mod$r[ , 1 ] )
    	SFs	<- cbind( SFs	, mod$s[ 1:2, ]	)
    	SFe	<- cbind( SFe	, mod$e[ 1:2, ]	)
    	if( boot == TRUE ){
    		mod_b	<-resf_boot( tau_sel = tau[ j ], q_sel = q[ j ], y = y,
    				     X = X, M = M, meig = meig, mod = mod, iter = iter )
    		SFb_boot_0  <- data.frame( t( mod_b$B ) )
    		SFs_boot_0  <- data.frame( t( mod_b$S ) )
    		SFb_boot2_00 <- t( apply( mod_b$B, 2, function( x ) quantile( x, probs = probs ) ) )
    		SFs_boot2_00 <- t( apply( mod_b$S, 2, function( x ) quantile( x, probs = probs ) ) )
    		SFb_p  <- 1 - abs( apply( mod_b$B, 2, function( x ) sum( x > 0 ) ) - iter / 2 ) / ( iter / 2 )
    		SFb_boot2_0  <- data.frame( mod$b[ , 1 ] , SFb_boot2_00, SFb_p )
    		SFs_boot2_0  <- data.frame( mod$s[ 1:2, ], SFs_boot2_00 )

    		rownames( SFb_boot_0 ) <- rownames( SFb_boot2_0 ) <- xname
    		rownames( SFs_boot_0 ) <- rownames( SFs_boot2_0 ) <- c( "shrink_sf_SE", "shrink_sf_alpha" )

    		names( SFb_boot_0 ) <- paste( "iter", 1:iter, sep = "" )
    		names( SFs_boot_0 ) <- paste( "iter", 1:iter, sep = "" )
    		SFb_boot[[ j ]]  <- SFb_boot_0
    		SFs_boot[[ j ]]  <- SFs_boot_0

    		names( SFb_boot2_0 ) <- c( "Estimates", "CI_lower", "CI_upper", "p_value" )
    		names( SFs_boot2_0 ) <- c( "Estimates", "CI_lower", "CI_upper" )
    		SFb_boot2[[ j ]]  <- SFb_boot2_0
    		SFs_boot2[[ j ]]  <- SFs_boot2_0

    		gc(); gc()
    	}
    	print( paste( "------- Complete: tau=", tau[ j ], " -------", sep = "" ) )
    }

    SFb		    <- data.frame( SFb )
    SFr		    <- data.frame( SFr )
    SFs		    <- data.frame( SFs )
    SFe		    <- data.frame( SFe )
    rownames( SFb ) <- xname
    rownames( SFr ) <- paste( "r", 1:ne, sep = "" )
    rownames( SFs ) <- c( "shrink_sf_SE", "shrink_sf_alpha" )
    rownames( SFe ) <- c( "resid_SE", "quasi_adjR2(cond)" )

    tau_name	        <- paste( "tau=", tau, sep= "" )
    names( SFb )      <- names( SFr ) <-names( SFs ) <- names( SFe ) <- tau_name
    if( boot == FALSE ){
      res   <- list( b = SFb, r = SFr, s = SFs, e = SFe, tau = tau )
    } else {
      tau_name2         <- paste( "--------------- tau=", tau, " ---------------", sep = "")
      names( SFb_boot ) <- names( SFb_boot2 ) <- tau_name2
      names( SFs_boot ) <- names( SFs_boot2 ) <- tau_name2
      res   <- list( b = SFb, r = SFr, s = SFs, e = SFe, B0 = SFb_boot,
                     S0 = SFs_boot, B = SFb_boot2, S = SFs_boot2, tau = tau )
      stopCluster( cl )
    }

return( res )
}


