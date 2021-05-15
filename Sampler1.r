library(expm)
##library(pracma)
CHAINS <- 10
STEPS <- 10
dt10 <- 0.0001
dt20 <- 0.0001
INTERVAL <- 1001
decaydt<-0.1
decayenergy<- 0.1
outbnd <- function(q) FALSE
clip1<-function(x,l,h){
    if(is.na(x)){
        -20
    }else{
        if(x<l){
            l
        }else{
            if(x>h)
                h
            else
                x
        }
    }
}

hmc <- function(U, dU, ddU, Dim, BURNIN, EPISODE, vanilla, switch,qinit){
    if(is.na(qinit)){
        qAll <- array(rnorm(CHAINS*Dim), dim=c(Dim,CHAINS))
    }else{
        qAll <- qinit
    }
    Utotal <- 0
    for(i in seq(1, CHAINS)) Utotal <- Utotal + U(qAll[,i])
    Htotal1 <- 2*Utotal
    Htotal2 <- 2*Utotal
    dt1 <- dt10
    dt2 <- dt20
    QS <- c()
    for(j in seq(1,EPISODE)){
        pAll <- array(rnorm(CHAINS*Dim), dim=c(Dim,CHAINS))
        KtotalNew <- 0
        for(i in seq(1, CHAINS)) {
            if(vanilla)
                KtotalNew <- KtotalNew + pAll[,i] %*% pAll[,i] / 2
            else
                KtotalNew <- KtotalNew + pAll[,i] %*% solve(ddU(qAll[,i]),pAll[,i]) /2
        }
        Utotal <- 0
        for(i in seq(1, CHAINS)) Utotal <- Utotal + U(qAll[,i])
        if(vanilla){
            Htotal <- Htotal1
            dt <- dt1
        }else{
            Htotal <- Htotal2
            dt <- dt2
        }
        Ktotal <- Htotal - Utotal
        pAll <- pAll * sqrt(abs(Ktotal/KtotalNew))[[1]]
        ES <- c()
        AS <- c()
        for(i in seq(1, CHAINS)) {
            bad <- FALSE
            p <- pAll[,i]
            q <- qAll[,i]
            UE <- c(U(q))
            q0 <- q
            for(s in seq(1, STEPS)) {
                p <- p - dt*dU(q)
                q1 <- q
                if(vanilla)
                    q <- q + dt*p
                else
                    q <- q + dt*solve(ddU(q), p)
                if(outbnd(q)){
                    q <- q1
                    bad <- TRUE
                }
                UE <- append(UE,U(q))
            }
            ES <- append(ES, UE)
            if(bad){
                alpha <- 0
            }else{
                #print(q0)
                #print(q)
                #print(U(q0))
                #print(U(q))
                alpha <- exp(clip1(U(q0)-U(q),-20,0))
            }
            AS <- append(AS, alpha)
            if(alpha < runif(1))
                q <- q0
            qAll[,i] <- q
            if(j>BURNIN)
                QS <- append(QS, q)
        }
        ES<-array(ES,dim=c(STEPS+1,CHAINS))
        s = c()
        S = c()
        for(i in seq(1, CHAINS)){
            s = append(s,order(ES[,i])[1])
            S = append(S,order(ES[,i])[STEPS+1])
        }
        s<-sort(unique(s))
        S<-sort(unique(S))
        if(j%%INTERVAL==0){
            print(j)
            #print(j,Ktotal,KtotalNew,mean(AS),sd(AS),Utotal,Htotal1,Htotal2,dt1,dt2,vanilla,s,S)
        }
        if(j<BURNIN){
            if(length(s)==2 & length(S)==2 & s[1]==1 & s[2]==STEPS+1 & S[1]==1 & S[2]==STEPS+1){
                dt = dt * (1+decaydt)
            }
            if(length(s)==1 & length(S)==1 & s[1]==1 & S[1]==STEPS+1){
                dt = dt / (1+decaydt)
            }
                if(mean(AS)>.9){
                    Htotal = (Htotal-Ktotal) * (1+decayenergy) + Utotal
                }else{
                    if(mean(AS)<.1)
                     Htotal = (Htotal-Ktotal) / (1+decayenergy) + Utotal
             
                }
        }
        if(vanilla){
            Htotal1 = Htotal
            dt1 = dt
        }else{
            Htotal2 = Htotal
            dt2 = dt
        }
        if(switch){
            vanilla=!vanilla
        }
    }
    array(QS,dim=c(Dim,CHAINS,EPISODE-BURNIN))
}


