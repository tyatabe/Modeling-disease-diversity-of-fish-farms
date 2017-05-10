setwd("C:/Users/tyatabe/OneDrive/Docs/PhD Epi/PhD project/Data/mortality")
# Getting score and disease data
score <- read.csv("# diseases.csv")


# Getting nw data for estimating centrality measures for each site
nw <- read.csv("Network_09_14.csv")
library(igraph)
el.14 <- graph.edgelist(as.matrix(nw[nw$Year== 2014,2:3]), directed=T)
el.13 <- graph.edgelist(as.matrix(nw[nw$Year==2013,2:3]), directed=T)
el.12 <- graph.edgelist(as.matrix(nw[nw$Year==2012,2:3]), directed=T)

# FHA-35: No movs in recorded during 2014 (in spite of reporting dis during survey) (only 2012: 4 movs in from 2 sites)
#FHA-36 fallowed since nov 2013 (1 mov in)...need to remove it (i.e. do nothing)
#FHA-69: no movs in on 2014, only movs on 2013(2 movs from 2 sites)
#FHA-47 no moves in on 2014 received. 2013: 7 in from 4 sites
# FHA 95 (Red mill), no recorded movs in: 0(?)

## adding an atribute for the edge to later have a simplified weighted nw
E(el.14)$link <- rep(1,141)
E(el.13)$link <- rep(1,131)
E(el.12)$link <- rep(1,172)

### Simplifying the graphs (i.e. removing multiple edges and loops)
el.14 <- simplify(el.14, edge.attr.comb="sum")
el.13 <- simplify(el.13, edge.attr.comb="sum")
el.12 <- simplify(el.12, edge.attr.comb="sum")

# Estimating centrality measures
indeg <-degree(el.14, v=V(el.14), mode="in")

indeg.w <- graph.strength (el.14, vids = V(el.14), mode = "in",
                           loops = F, weights = E(el.14)$link)

indeg.w13 <- graph.strength (el.13, vids = V(el.13), mode = "in",
                           loops = F, weights = E(el.13)$link)
indeg.w12 <- graph.strength (el.12, vids = V(el.12), mode = "in",
                           loops = F, weights = E(el.12)$link)


# Data frame of site's windeg for each year
Windeg <- unname(c(indeg.w, indeg.w13, indeg.w12))
FHA <- c(attr(indeg.w, "names"), attr(indeg.w13, "names"), 
                   attr(indeg.w12, "names"))
Year <- c(rep(2014, 65), rep(2013, 62), rep(2012, 77))


Suppliers <- as.data.frame(cbind(Windeg, FHA, Year))
write.csv(Suppliers, file="Suppliers.csv")


# Adding sites' degree with data from 2013 manually
deg.m <- c(2, 2, 4, 0)
windeg.m <- c(4, 2, 7, 0)

#deg.m <- c(0, 0)
#windeg.m <- c(0, 0)


fha.m <- c("FHA-35", "FHA-69", "FHA-47", "FHA-95")

# fha.m <- c("FHA-35", "FHA-95")

### Creating a data frame with farm centrality measures + data for sites not moving in 2014
Centrality<-as.data.frame(c(indeg, deg.m))
colnames(Centrality)<- "indeg"
Centrality$windeg <- c(indeg.w, windeg.m)
Centrality$ID <- as.factor(c(attr(indeg, "names"), fha.m))




# Adding indegree to score dataset
dat <- merge(score, Centrality, by.x= "FHA", by.y="ID")
# write .csv of data
write.csv(dat, file="dat.csv")

# Done with igraph, detach it
detach("package:igraph", unload=TRUE)

# Subsetting only to SW sites
d <- dat[dat$Type=="SW",]
d2 <- dat[dat$Type!="SW",]
d3 <- dat
d3$Type <- as.factor(ifelse(d3$Type=="SW", "SW", ifelse(d3$Type=="FW-trout", "Trout","FW") ))


# Descriptive stats of score and indegree
library(doBy)
cuart1 <- function(x){quantile(x, probs=0.25)}
cuart2 <- function(x){quantile(x, probs=0.50)}
cuart3 <- function(x){quantile(x, probs=0.75)}


summaryBy(Score_surv100 ~ Type , FUN=c(cuart1, cuart2, cuart3), data=score)
summaryBy(indeg  ~ Type, FUN=c(cuart1, cuart2, cuart3, mean), data=dat)
summaryBy(No_outbreaks  ~ Type, FUN=c(cuart1, cuart2, cuart3, max), data=dat)


### Fitting bayesian Poisson models
library(rethinking)

# Remove variables with missing values (otherwiese stan models crash) and variables not used
d <- d[,c("Score_surv100", "Score_55", "Manager", "No_inf", "No_outbreaks", "indeg", "windeg")]
d2 <- d2[,c("Score_surv100", "Score_55", "Manager", "No_inf", "No_outbreaks", "indeg", "windeg")]
d3 <- d3[,c("Type","Score_surv100", "Score_55", "Manager", "No_inf", "No_outbreaks", "indeg", "windeg")]
d3$Score_surv100 <-ifelse(d3$Type=="FW", d3$Score_surv100/max(d3[d3$Type=="FW","Score_surv100"]),
               ifelse(d3$Type=="Trout", d3$Score_surv100/max(d3[d3$Type=="Trout","Score_surv100"]),
                      d3$Score_surv100/max(d3[d3$Type=="SW","Score_surv100"]))) 



# Centering
# construct scaled (by 2sd) score and centered (minus mean) categorized indeg (1 if >1 indeg)
d$windeg_c <- ifelse(d$indeg>1, 1,0)
d$windeg_c <- d$windeg_c-mean(d$windeg_c)
d$Score_surv100_c <- scale(d$Score_surv100, scale=sd(d$Score_surv100)*2)
d$inter <- d$Score_surv100_c*d$windeg_c
d$Score_55_c <- d$Score_55 - mean(d$Score_55)

d2$windeg_c <- ifelse(d2$indeg>1, 1,0)
d2$windeg_c <- d2$windeg_c-mean(d2$windeg_c)
d2$Score_surv100_c <- scale(d2$Score_surv100, scale=sd(d2$Score_surv100)*2)
d2$inter <- d2$Score_surv100_c*d2$windeg_c
d2$Score_55_c <- d2$Score_55 - mean(d2$Score_55)

d3$windeg_c <- ifelse(d3$indeg>1, 1,0)
d3$windeg_c <- d3$windeg_c-mean(d3$windeg_c)
d3$Score_surv100_c <- scale(d3$Score_surv100, scale=sd(d3$Score_surv100)*2)
d3$inter <- d3$Score_surv100_c*d3$windeg_c
d3$Score_55_c <- d3$Score_55 - mean(d3$Score_55)


## re-estimate
# Full model
m5.1 <- map2stan(
alist(No_outbreaks ~ dpois( lambda ),
    log(lambda) <- a + bs*Score_surv100_c+ bi*windeg_c + bsi*Score_surv100_c*windeg_c,
    a ~ dnorm(0,100),
    c(bs,bi,bsi) ~ dnorm(0,1)
  ),
  data=d, iter=3000 , warmup=1000 , chains=4 ,
    control = list(adapt_delta = 0.80))

precis(m5.1, prob=0.95, digits=3)
pairs(m5.1)

# No interaction
m5.2 <- map2stan(
alist(No_outbreaks ~ dpois( lambda ),
    log(lambda) <- a + bs*Score_surv100_c+ bi*windeg_c ,
    a ~ dnorm(0,100),
    c(bs,bi) ~ dnorm(0,1)
  ),
  data=d, iter=20000 , warmup=1000 , chains=1 ,
    control = list(adapt_delta = 0.80), start=list(bs=0))

precis(m5.2, prob=0.95, digits=3)
pairs(m5.2)


# No windeg
m5.3 <- map2stan(
alist(No_outbreaks ~ dpois( lambda ),
    log(lambda) <- a + bs*Score_surv100_c,
    a ~ dnorm(0,100),
    c(bs) ~ dnorm(0,1)
  ),
  data=d, iter=20000 , warmup=1000 , chains=1 )

precis(m5.3, prob=0.95, digits=3)
pairs(m5.3)

# No biosec
m5.4 <- map2stan(
alist(No_outbreaks ~ dpois( lambda ),
    log(lambda) <- a + bi*windeg_c ,
    a ~ dnorm(0,100),
    c(bi) ~ dnorm(0,1)
  ),
  data=d, iter=20000 , warmup=1000 , chains=1 ,
control = list(adapt_delta = 0.80))

precis(m5.4, prob=0.95, digits=3)
pairs(m5.4)

# Null model
m5.5 <- map2stan(
alist(No_outbreaks ~ dpois( lambda ),
    log(lambda) <- a ,
    a ~ dnorm(0,100)
  ),
  data=d, iter=20000 , warmup=1000 , chains=1 ,
control = list(adapt_delta = 0.80))

precis(m5.5, prob=0.95, digits=3)



# Comparing models

diseases.compare.stan <- compare(m5.1,m5.2,m5.3,m5.4,m5.5,n=1e4)
plot(diseases.compare.stan)

# Comparing parameter values accross models
coeftab(m5.1 , m5.2 , m5.3, m5.4, m5.5)
plot(coeftab(m5.1 , m5.2 , m5.3, m5.4, m5.5))



## Comparing high and low contact for fixed high and low biosec
post5.1 <- extract.samples(m5.1)

# Low biosec
dis_high_c <- exp(post5.1$a + post5.1$bi + (post5.1$bsi + post5.1$bs)*-2.4006)
dis_low_c <- exp(post5.1$a + post5.1$bs*-2.4006)

# Difference in number of diseases
dis_diff <- dis_high_c - dis_low_c
sum(dis_diff > 0)/length(dis_diff)
plot(density(dis_diff), xlab="Difference high vs low contact|low biosec", main="",
     xlim=c(-4, 4))
shade(density(dis_diff), HPDI(dis_diff))
abline(v=median(dis_diff), lty=4)

# High biosec
dis_high_c <- exp(post5.1$a + post5.1$bi + (post5.1$bsi + post5.1$bs)*3.0432)
dis_low_c <- exp(post5.1$a + post5.1$bs*3.0432)

# Difference in number of diseases
dis_diff <- dis_high_c - dis_low_c
sum(dis_diff > 0)/length(dis_diff)
plot(density(dis_diff), xlab="Difference high vs low contact|high biosec", main="",
     xlim=c(-4, 4))
shade(density(dis_diff), HPDI(dis_diff))
abline(v=median(dis_diff), lty=4)

## Comparing high and low biosec for fixed high and low contact rate

# High contact rate
dis_high_bs <- exp(post5.1$a + post5.1$bi + (post5.1$bsi + post5.1$bs)*3.0432)
dis_low_bs <- exp(post5.1$a + post5.1$bi + (post5.1$bsi + post5.1$bs)*-2.4006)

# Difference in number of diseases
dis_diff <- dis_low_bs - dis_high_bs
sum(dis_diff > 0)/length(dis_diff)
plot(density(dis_diff), xlab="Difference low vs high biosecurity|high contact", main="",
     xlim=c(-4, 4))
shade(density(dis_diff), HPDI(dis_diff))
abline(v=median(dis_diff), lty=4)

# Low contact rate
dis_high_bs <- exp(post5.1$a + post5.1$bs* 3.0432)
dis_low_bs <- exp(post5.1$a + post5.1$bs*-2.4006)

# Difference in number of diseases
dis_diff <- dis_low_bs - dis_high_bs
sum(dis_diff > 0)/length(dis_diff)
plot(density(dis_diff), xlab="Difference low vs high biosecurity|low contact", main="",
     xlim=c(-4, 4))
shade(density(dis_diff), HPDI(dis_diff))
abline(v=median(dis_diff), lty=4)



# Doing it with model averaging

# High contact vs low contact given low biosec

Score_surv100_c <- -2.4006
windeg_c <- 1
d.pred_hc_lbs <- data.frame(Score_surv100_c, windeg_c)

lambda.pred.hc_lbs <- ensemble( m5.1 , m5.2 , m5.3, m5.4, m5.5 , data=d.pred_hc_lbs )

Score_surv100_c <- -2.4006
windeg_c <- 0
d.pred_lc_lbs <- data.frame(Score_surv100_c, windeg_c)

lambda.pred.lc_lbs <- ensemble( m5.1 , m5.2 , m5.3, m5.4, m5.5 , data=d.pred_lc_lbs )


# Difference in number of diseases
dis_diff <- lambda.pred.hc_lbs$link - lambda.pred.lc_lbs$link
sum(dis_diff > 0)/length(dis_diff)
plot(density(dis_diff), xlab="Difference high vs low contact|low biosec", main="",
     xlim=c(-4, 4), ylim=c(0, 1.1))
shade(density(dis_diff), PI(dis_diff, prob=0.90))
abline(v=median(dis_diff), lty=2)


# Top model
dis_high_c <- exp(post5.1$a + post5.1$bi + (post5.1$bsi + post5.1$bs)*-2.4006)
dis_low_c <- exp(post5.1$a + post5.1$bs*-2.4006)

# Difference in number of diseases
dis_diff <- dis_high_c - dis_low_c
sum(dis_diff > 0)/length(dis_diff)
lines(density(dis_diff), lty=5)
shade(density(dis_diff), PI(dis_diff, prob=0.9), col=col.alpha(rangi2,0.2))
abline(v=median(dis_diff), lty=5)


# High contact vs low contact given high biosec

Score_surv100_c <- 3.0432
windeg_c <- 1
d.pred_hc_hbs <- data.frame(Score_surv100_c, windeg_c)

lambda.pred.hc_hbs <- ensemble( m5.1 , m5.2 , m5.3, m5.4, m5.5 , data=d.pred_hc_hbs )

Score_surv100_c <- 3.0432
windeg_c <- 0
d.pred_lc_hbs <- data.frame(Score_surv100_c, windeg_c)

lambda.pred.lc_hbs <- ensemble( m5.1 , m5.2 , m5.3, m5.4, m5.5 , data=d.pred_lc_hbs )


# Difference in number of diseases
dis_diff <- lambda.pred.hc_hbs$link - lambda.pred.lc_hbs$link
sum(dis_diff > 0)/length(dis_diff)
plot(density(dis_diff, adjust=3), xlab="Difference high vs low contact|low biosec", main="",
     xlim=c(-4, 4), ylim=c(0, 0.7))
shade(density(dis_diff, adjust=3), PI(dis_diff, prob=0.9))
abline(v=median(dis_diff), lty=2)


# Top model
dis_high_c <- exp(post5.1$a + post5.1$bi + (post5.1$bsi + post5.1$bs)*3.0432)
dis_low_c <- exp(post5.1$a + post5.1$bs*3.0432)

# Difference in number of diseases
dis_diff <- dis_high_c - dis_low_c
sum(dis_diff > 0)/length(dis_diff)
lines(density(dis_diff), lty=5)
shade(density(dis_diff), PI(dis_diff, prob=0.9), col=col.alpha(rangi2,0.2))
abline(v=median(dis_diff), lty=5)




# High biosec vs low biosec given high contact

Score_surv100_c <- 3.0432
windeg_c <- 1
d.pred_hbs_hc <- data.frame(Score_surv100_c, windeg_c)

lambda.pred.hbs_hc <- ensemble( m5.1 , m5.2 , m5.3, m5.4, m5.5 , data=d.pred_hbs_hc)

Score_surv100_c <- -2.4006
windeg_c <- 1
d.pred_lbs_hc <- data.frame(Score_surv100_c, windeg_c)

lambda.pred.lbs_hc <- ensemble( m5.1 , m5.2 , m5.3, m5.4, m5.5 , data=d.pred_lbs_hc)


# Difference in number of diseases
dis_diff <- lambda.pred.lbs_hc$link - lambda.pred.hbs_hc$link
sum(dis_diff > 0)/length(dis_diff)
plot(density(dis_diff, adjust=2), xlab="Difference low vs high biosecurity|high contact", main="",
     xlim=c(-4, 4), ylim=c(0, 1.0))
shade(density(dis_diff, adjust=2), PI(dis_diff, prob=0.9))
abline(v=median(dis_diff), lty=2)

dis_high_bs <- exp(post5.1$a + post5.1$bi + (post5.1$bsi + post5.1$bs)*3.0432)
dis_low_bs <- exp(post5.1$a + post5.1$bi + (post5.1$bsi + post5.1$bs)*-2.4006)

# Difference in number of diseases
dis_diff <- dis_low_bs - dis_high_bs
sum(dis_diff > 0)/length(dis_diff)
lines(density(dis_diff), lty=5)
shade(density(dis_diff), PI(dis_diff, prob=0.9), col=col.alpha(rangi2,0.2))
abline(v=median(dis_diff), lty=5)




# High biosec vs low biosec given low contact

Score_surv100_c <- 3.0432
windeg_c <- 0
d.pred_hbs_lc <- data.frame(Score_surv100_c, windeg_c)

lambda.pred.hbs_lc <- ensemble( m5.1 , m5.2 , m5.3, m5.4, m5.5 , data=d.pred_hbs_lc)

Score_surv100_c <- -2.4006
windeg_c <- 0
d.pred_lbs_lc <- data.frame(Score_surv100_c, windeg_c)

lambda.pred.lbs_lc <- ensemble( m5.1 , m5.2 , m5.3, m5.4, m5.5 , data=d.pred_lbs_lc)


# Difference in number of diseases
dis_diff <- lambda.pred.lbs_lc$link - lambda.pred.hbs_lc$link
sum(dis_diff > 0)/length(dis_diff)
plot(density(dis_diff, adjust=1), xlab="Difference low vs high biosecurity|low contact", main="",
     xlim=c(-4, 4), ylim=c(0, 1.2))
shade(density(dis_diff, adjust=1), PI(dis_diff, prob=0.9))
abline(v=median(dis_diff), lty=2)


dis_high_bs <- exp(post5.1$a + post5.1$bs* 3.0432)
dis_low_bs <- exp(post5.1$a + post5.1$bs*-2.4006)



# Difference in number of diseases
dis_diff <- dis_low_bs - dis_high_bs
sum(dis_diff > 0)/length(dis_diff)
lines(density(dis_diff), lty=5)
shade(density(dis_diff), PI(dis_diff, prob=0.9), col=col.alpha(rangi2,0.2))
abline(v=median(dis_diff), lty=5)



# Evaluating the model in predicting (retrodicting) oberseved data 
sim.obs <- sim(m5.1, data=d)
sim.obs.90 <- apply(sim.obs, 2, PI, prob=0.9)
sim.obs.ens <- ensemble( m5.1 , m5.2 , m5.3, m5.4, m5.5 , data=d )
sim.obs.ens.90 <- apply(sim.obs.ens$sim, 2, PI, prob=0.9)


check <- postcheck(m5.1, prob=0.9)
points(seq(1:17), apply(sim.obs, 2, median), pch=16, col="red")

sim.rate <- link(m5.1, data=d)
points(seq(1:17), rpois(17, apply(sim.rate, 2, median)), pch=16, col="red")


# Do counterfactual predictions

# sequence of biosecurity scores and biosecs to compute over
biosec.seq <- seq( from=-5.5 , to=5.5 , length.out=80 )

# compute trend for high contact rate site
d.pred_hc <- data.frame(
  Score_surv100_c = biosec.seq,
  windeg_c = 1
)

lambda.pred.h <- ensemble( m5.1 , m5.2 , m5.3, m5.4, m5.5 , data=d.pred_hc )
lambda.med <- apply( lambda.pred.h$link , 2 , median )
lambda.PI <- apply( lambda.pred.h$link , 2 , PI, prob=0.66 )
lambda.PI.90 <- apply( lambda.pred.h$link , 2 , PI, prob=0.90 )

# compute trend for low contact rate site
d.pred_lc <- data.frame(
  Score_surv100_c = biosec.seq,
  windeg_c = 0
)

lambda.pred.l <- ensemble( m5.1 , m5.2 , m5.3, m5.4, m5.5 , data=d.pred_lc )
lambda.med.l <- apply( lambda.pred.l$link , 2 , median )
lambda.PI.l <- apply( lambda.pred.l$link , 2 , PI, prob=0.66 )
lambda.PI.l.90 <- apply( lambda.pred.l$link , 2 , PI, prob=0.90 )


# Plotting
pch <- ifelse( d$indeg > 1 , 16 , 1 )
set.seed(0)
plot( jitter(d$Score_surv100_c, factor=15) , jitter(d$No_outbreaks, factor=0.0) , col=rangi2 , pch=pch ,
      xlab="Biosecurity score" , ylab="No of different diseases", ylim=c(0,6), xlim=c(-5.5, 5.5) )

# plot predicted rates
lines( biosec.seq , lambda.med , col=rangi2 )
shade( lambda.PI.90 , biosec.seq , col=col.alpha("pink",0.2) )
shade( lambda.PI , biosec.seq , col=col.alpha(rangi2,0.2) )



lines( biosec.seq , lambda.med.l , lty=2 )
shade( lambda.PI.l.90 , biosec.seq , col=col.alpha("green4",0.2))
shade( lambda.PI.l , biosec.seq , col=col.alpha("black",0.1) )


# Now without model averaging/ensemble
# Model 5.1
lambda.pred.5.1.h <- link(m5.1, data=d.pred_hc)
lambda.med.5.1 <- apply( lambda.pred.5.1.h, 2 , median )
lambda.PI5.1 <- apply( lambda.pred.5.1.h, 2 , PI, prob=0.66 )
lambda.PI5.1.90 <- apply( lambda.pred.5.1.h, 2 , PI, prob=0.90 )


lambda.pred.5.1.l <- link(m5.1, data=d.pred_lc)
lambda.med.5.1.l <- apply( lambda.pred.5.1.l, 2 , median )
lambda.PI5.1.l <- apply( lambda.pred.5.1.l, 2 , PI, prob=0.66 )
lambda.PI5.1.l.90 <- apply( lambda.pred.5.1.l, 2 , PI, prob=0.90 )


# Plotting
set.seed(0)
plot( jitter(d$Score_surv100_c, factor=15) , jitter(d$No_outbreaks, factor=0.0) , col=rangi2 , pch=pch ,
      xlab="Biosecurity score" , ylab="No of different diseases", ylim=c(0,6), xlim=c(-5.5,5.5) )

lines( biosec.seq , lambda.med.5.1 , col=rangi2 )
shade( lambda.PI5.1.90 , biosec.seq , col=col.alpha("pink",0.2) )
shade( lambda.PI5.1 , biosec.seq , col=col.alpha(rangi2,0.2) )

lines( biosec.seq , lambda.med.5.1.l , lty=2 )
shade( lambda.PI5.1.l.90 , biosec.seq , col=col.alpha("green",0.1) )
shade( lambda.PI5.1.l , biosec.seq , col=col.alpha("black",0.1) )



# Plotting predictions

lambda.sim.5.1.h <- sim(m5.1, data=d.pred_hc)
lambda.sim.med.5.1 <- apply( lambda.sim.5.1.h, 2 , median )
lambda.sim.PI5.1 <- apply( lambda.sim.5.1.h, 2 , PI, prob=0.66 )
lambda.sim.PI5.1.90 <- apply( lambda.sim.5.1.h, 2 , PI, prob=0.90 )


lambda.sim.5.1.l <- sim(m5.1, data=d.pred_lc)
lambda.sim.med.5.1.l <- apply( lambda.sim.5.1.l, 2 , median )
lambda.sim.PI5.1.l <- apply( lambda.sim.5.1.l, 2 , PI, prob=0.66 )
lambda.sim.PI5.1.l.90 <- apply( lambda.sim.5.1.l, 2 , PI, prob=0.90 )




set.seed(0)
plot( jitter(d$Score_surv100_c, factor=15) , jitter(d$No_outbreaks, factor=0.0) , col=rangi2 , pch=pch ,
      xlab="Biosecurity score" , ylab="No of different diseases", ylim=c(0,6), xlim=c(-5.5,5.5) )

lines( biosec.seq , lambda.sim.med.5.1 , col=rangi2 )
shade( lambda.sim.PI5.1.90 , biosec.seq , col=col.alpha("pink",0.2) )
shade( lambda.sim.PI5.1 , biosec.seq , col=col.alpha(rangi2,0.2) )

lines( biosec.seq , lambda.sim.med.5.1.l , lty=2 )
shade( lambda.sim.PI5.1.l.90 , biosec.seq , col=col.alpha("green",0.1) )
shade( lambda.sim.PI5.1.l , biosec.seq , col=col.alpha("black",0.1) )





# Model for number of infectious diseases 


# Full model
m6.1 <- map2stan(
  alist(No_inf ~ dpois( lambda ),
        log(lambda) <- a + bs*Score_surv100_c + bi*windeg_c + bsi*Score_surv100_c*windeg_c,
        a ~ dnorm(0,100),
        c(bs,bi,bsi) ~ dnorm(0,1)
  ),
  data=d, iter=20000 , warmup=1000 , chains=1 ,
  control = list(adapt_delta = 0.80))

precis(m6.1, prob=0.95, digits=3)
pairs(m6.1)

# No interaction
m6.2 <- map2stan(
  alist(No_inf ~ dpois( lambda ),
        log(lambda) <- a + bs*Score_surv100_c + bi*windeg_c ,
        a ~ dnorm(0,100),
        c(bs,bi) ~ dnorm(0,1)
  ),
  data=d, iter=20000 , warmup=1000 , chains=1 ,
  control = list(adapt_delta = 0.80))

precis(m6.2, prob=0.95, digits=3)
pairs(m6.2)


# No windeg
m6.3 <- map2stan(
  alist(No_inf ~ dpois( lambda ),
        log(lambda) <- a + bs*Score_surv100_c,
        a ~ dnorm(0,100),
        c(bs) ~ dnorm(0,1)
  ),
  data=d, iter=20000 , warmup=1000 , chains=1 )

precis(m6.3, prob=0.95, digits=3)
pairs(m6.3)

# No biosec
m6.4 <- map2stan(
  alist(No_inf ~ dpois( lambda ),
        log(lambda) <- a + bi*windeg_c ,
        a ~ dnorm(0,100),
        c(bi) ~ dnorm(0,1)
  ),
  data=d, iter=20000 , warmup=1000 , chains=1 ,
  control = list(adapt_delta = 0.80))

precis(m6.4, prob=0.95, digits=3)
pairs(m6.4)

# Null model
m6.5 <- map2stan(
  alist(No_inf ~ dpois( lambda ),
        log(lambda) <- a ,
        a ~ dnorm(0,100)
  ),
  data=d, iter=20000 , warmup=1000 , chains=1 ,
  control = list(adapt_delta = 0.80))

precis(m6.5, prob=0.95, digits=3)
plot(m6.5)


# Comparing models

diseases.compare.stan_dis <- compare(m6.1,m6.2,m6.3,m6.4,m6.5,n=1e4)
plot(diseases.compare.stan_dis)


coeftab(m6.1,m6.2,m6.3,m6.4,m6.5)
plot(coeftab(m6.1,m6.2,m6.3,m6.4,m6.5))


# Check implied predictions for high contact rate
d.pred_hc 
d.pred_lc


lambda.pred.h <- ensemble( m6.1 , m6.2 , m6.3, m6.4, m6.5 , data=d.pred_hc )
lambda.med.h <- apply( lambda.pred.h$link , 2 , median )
lambda.PI.h <- apply( lambda.pred.h$link , 2 , HPDI )

# Plotting trend for high contact sites
set.seed(0)
plot( jitter(d$Score_surv100_c, factor=15) , jitter(d$No_inf, factor=0.0) , col=rangi2 , pch=pch ,
      xlab="Biosecurity score" , ylab="No of different diseases", ylim=c(0,3.5), xlim=c(-5.5,5.5) )
lines( biosec.seq , lambda.med.h , col=rangi2 )
shade( lambda.PI.h , biosec.seq , col=col.alpha(rangi2,0.2) )

# compute trend for low contact sites


lambda.pred.l <- ensemble( m6.1 , m6.2 , m5.3, m6.4, m5.5 , data=d.pred_lc )
lambda.med.l <- apply( lambda.pred.l$link , 2 , median )
lambda.PI.l <- apply( lambda.pred.l$link , 2 , HPDI )

# Plotting trend for high contact sites
lines( biosec.seq , lambda.med.l , lty=2 )
shade( lambda.PI.l , biosec.seq , col=col.alpha("black",0.1) )


# Checking for model 6.3
# 6.3 high contact

lambda.pred.6.3_h <- link(m6.3, data=d.pred_hc)
lambda.med.6.3.h <- apply( lambda.pred.6.3_h, 2 , median )
lambda.PI.6.3.h <- apply( lambda.pred.6.3_h, 2 , HPDI )

set.seed(0)
plot( jitter(d$Score_surv100_c, factor=15) , jitter(d$No_inf, factor=0.5) , col=rangi2 , pch=16 ,
      xlab="Biosecurity score" , ylab="No of different diseases", ylim=c(0,3.5), xlim=c(-5.5,5.5) )
lines( biosec.seq , lambda.med.6.3.h , col=rangi2 )
shade( lambda.PI.6.3.h , biosec.seq , col=col.alpha(rangi2,0.2) )



# Differences in number of diseases
post6.3 <- extract.samples(m6.3)

dis_high_bs <- exp(post6.3$a + post6.3$bs*3.3093)
dis_low_bs <- exp(post6.3$a + post6.3$bs*-2.3218)

# Difference in number of diseases
dis_diff <- dis_low_bs - dis_high_bs
sum(dis_diff > 0)/length(dis_diff)
plot(density(dis_diff), xlab="Difference low vs high biosecurity", main="")
shade(density(dis_diff), HPDI(dis_diff))
abline(v=median(dis_diff), lty=4)


# Fitting model with score based on 55 variables only


# Full model
m7.1 <- map2stan(
  alist(No_outbreaks ~ dpois( lambda ),
        log(lambda) <- a + bs*Score_55_c+ bi*windeg_c + bsi*Score_55_c*windeg_c,
        a ~ dnorm(0,100),
        c(bs,bi,bsi) ~ dnorm(0,1)
  ),
  data=d, iter=20000 , warmup=1000 , chains=1 ,
  control = list(adapt_delta = 0.80))

precis(m7.1, prob=0.95, digits=3)
pairs(m7.1)

# No interaction
m7.2 <- map2stan(
  alist(No_outbreaks ~ dpois( lambda ),
        log(lambda) <- a + bs*Score_55_c+ bi*windeg_c ,
        a ~ dnorm(0,100),
        c(bs,bi) ~ dnorm(0,1)
  ),
  data=d, iter=20000 , warmup=1000 , chains=1 ,
  control = list(adapt_delta = 0.80))

precis(m7.2, prob=0.95, digits=3)
pairs(m7.2)


# No windeg
m7.3 <- map2stan(
  alist(No_outbreaks ~ dpois( lambda ),
        log(lambda) <- a + bs*Score_55_c,
        a ~ dnorm(0,100),
        c(bs) ~ dnorm(0,1)
  ),
  data=d, iter=20000 , warmup=1000 , chains=1 )

precis(m7.3, prob=0.95, digits=3)
pairs(m7.3)

# No biosec
m7.4 <- map2stan(
  alist(No_outbreaks ~ dpois( lambda ),
        log(lambda) <- a + bi*windeg_c ,
        a ~ dnorm(0,100),
        c(bi) ~ dnorm(0,1)
  ),
  data=d, iter=20000 , warmup=1000 , chains=1 ,
  control = list(adapt_delta = 0.80))

precis(m7.4, prob=0.95, digits=3)
pairs(m7.4)

# Null model
m7.5 <- map2stan(
  alist(No_outbreaks ~ dpois( lambda ),
        log(lambda) <- a ,
        a ~ dnorm(0,100)
  ),
  data=d, iter=20000 , warmup=1000 , chains=1 ,
  control = list(adapt_delta = 0.80))

precis(m7.5, prob=0.95, digits=3)



# Comparing models

diseases.compare.stan_55 <- compare(m7.1,m7.2,m7.3,m7.4,m7.5,n=1e4)
plot(diseases.compare.stan_55)

# Comparing models based on different scores (full vs reduced)
compare(m5.1,m5.2,m5.3,m5.4,m5.5, m7.1,m7.2,m7.3,m7.4,m7.5,n=1e4)
plot(compare(m5.1,m5.2,m5.3,m5.4,m5.5, m7.1,m7.2,m7.3,m7.4,m7.5,n=1e4))


## Comparing high and low contact for fixed high and low biosec
post5.1 <- extract.samples(m5.1)

# Low biosec
dis_high_c <- exp(post5.1$a + post5.1$bi + (post5.1$bsi + post5.1$bs)*-2.3218)
dis_low_c <- exp(post5.1$a + post5.1$bs*-2.3218)

# Difference in number of diseases
dis_diff <- dis_high_c - dis_low_c
sum(dis_diff > 0)/length(dis_diff)
plot(density(dis_diff), xlab="Difference high vs low contact|low biosec", main="",
     xlim=c(-4, 4))
shade(density(dis_diff), HPDI(dis_diff))
abline(v=median(dis_diff), lty=4)

# High biosec
dis_high_c <- exp(post5.1$a + post5.1$bi + (post5.1$bsi + post5.1$bs)*3.3093)
dis_low_c <- exp(post5.1$a + post5.1$bs*3.3093)

# Difference in number of diseases
dis_diff <- dis_high_c - dis_low_c
sum(dis_diff > 0)/length(dis_diff)
plot(density(dis_diff), xlab="Difference high vs low contact|high biosec", main="",
     xlim=c(-4, 4))
shade(density(dis_diff), HPDI(dis_diff))
abline(v=median(dis_diff), lty=4)

## Comparing high and low biosec for fixed high and low contact rate

# High contact rate
dis_high_bs <- exp(post5.1$a + post5.1$bi + (post5.1$bsi + post5.1$bs)*3.3093)
dis_low_bs <- exp(post5.1$a + post5.1$bi + (post5.1$bsi + post5.1$bs)*-2.3218)

# Difference in number of diseases
dis_diff <- dis_low_bs - dis_high_bs
sum(dis_diff > 0)/length(dis_diff)
plot(density(dis_diff), xlab="Difference low vs high biosecurity|high contact", main="",
     xlim=c(-4, 4))
shade(density(dis_diff), HPDI(dis_diff))
abline(v=median(dis_diff), lty=4)

# Low contact rate
dis_high_bs <- exp(post5.1$a + post5.1$bs*3.3093)
dis_low_bs <- exp(post5.1$a + post5.1$bs*-2.3218)

# Difference in number of diseases
dis_diff <- dis_low_bs - dis_high_bs
sum(dis_diff > 0)/length(dis_diff)
plot(density(dis_diff), xlab="Difference low vs high biosecurity|low contact", main="",
     xlim=c(-4, 4))
shade(density(dis_diff), HPDI(dis_diff))
abline(v=median(dis_diff), lty=4)


# Do counterfactual predictions
# make plot of raw data to begin
pch <- ifelse( d$indeg > 1 , 16 , 1 )
set.seed(0)
plot( jitter(d$Score_surv100_c, factor=3) , jitter(d$No_outbreaks, factor=0.5) , col=rangi2 , pch=pch ,
      xlab="Biosecurity score" , ylab="No of different diseases", ylim=c(0,6) )

# sequence of biosecurity scores and biosecs to compute over
biosec.seq <- seq( from=-5.5 , to=5 , length.out=80 )

# compute trend for high contact rate site
d.pred_hc <- data.frame(
  Score_surv100_c = biosec.seq,
  windeg_c = 1
)

lambda.pred.h <- ensemble( m5.1 , m5.2 , m5.3, m5.4, m5.5 , data=d.pred_hc )
lambda.med <- apply( lambda.pred.h$link , 2 , median )



lambda.PI <- apply( lambda.pred.h$link , 2 , HPDI )

# plot predicted trend
lines( biosec.seq , lambda.med , col=rangi2 )
shade( lambda.PI , biosec.seq , col=col.alpha(rangi2,0.2) )


# compute trend for low contact rate site
d.pred_lc <- data.frame(
  Score_surv100_c = biosec.seq,
  windeg_c = 0
)

lambda.pred.h <- ensemble( m5.1 , m5.2 , m5.3, m5.4, m5.5 , data=d.pred_lc )
lambda.med <- apply( lambda.pred.h$link , 2 , median )
lambda.PI <- apply( lambda.pred.h$link , 2 , HPDI )

lines( biosec.seq , lambda.med , lty=2 )
shade( lambda.PI , biosec.seq , col=col.alpha("black",0.1) )



# Now without model averaging/ensemble
# Model 5.1
lambda.pred.5.1.h <- link(m5.1, data=d.pred_hc)
lambda.med.5.1 <- apply( lambda.pred.5.1.h, 2 , median )
lambda.PI5.1 <- apply( lambda.pred.5.1.h, 2 , HPDI )

lambda.pred.5.1.l <- link(m5.1, data=d.pred_lc)
lambda.med.5.1.l <- apply( lambda.pred.5.1.l, 2 , median )
lambda.PI5.1.l <- apply( lambda.pred.5.1.l, 2 , HPDI )


# Plotting
set.seed(0)
plot( jitter(d$Score_surv100_c, factor=3) , jitter(d$No_outbreaks, factor=0.5) , col=rangi2 , pch=pch ,
      xlab="Biosecurity score" , ylab="No of different diseases", ylim=c(0,6.5) )

lines( biosec.seq , lambda.med.5.1 , col=rangi2 )
shade( lambda.PI5.1 , biosec.seq , col=col.alpha(rangi2,0.2) )

lines( biosec.seq , lambda.med.5.1.l , lty=2 )
shade( lambda.PI5.1.l , biosec.seq , col=col.alpha("black",0.1) )






# Need to do a PCA for score variables of SW sites and then regress on the most important
# components(?)
# Read in data
setwd("C:/Users/tyatabe/OneDrive/Docs/PhD Epi/PhD project/Data/Biosec_Survey")
d2 <- read.csv("PCA_SW_csv.csv")
# Centering(standardizing?) data

# PCA
pca <- prcomp(d2[,-1], scale=T)

# Plotting the scree plot and cum variance plor
pve =100* pca$sdev ^2/ sum (pca$sdev ^2)

plot(pve , type ="o", ylab =" PVE ", xlab =" Principal Component ",
       col =" blue ")
plot( cumsum (pve ), type ="o", ylab =" Cumulative PVE ", xlab ="
Principal Component ", col =" brown3 ")

# Getting the first 5 components' loadings
loads <- pca$rotation[,1:5]
# Creating a csv file to check on excel
write.csv(loads, "loads_SW.csv")

# Plotting the components
plot(pca$x[,1:2], pch=16)
abline(h=0, v=0, col="grey", lty=2)
text(pca$x[,1:2], labels=score[1:18,11], pos=1)

# Plotting the components
plot(pca$x[,3:4], pch=16)
abline(h=0, v=0, col="grey", lty=2)
text(pca$x[,3:4], labels=score[1:18,11], pos=1)

# Visually checking if components are related to # dis
comps <- data.frame(score[1:18,c(2,11)], pca$x[,1:5])
write.csv(comps, "comps_SW.csv")


# Fitting a model with score based on 49 variables only
setwd("C:/Users/tyatabe/OneDrive/Docs/PhD Epi/PhD project/Data/mortality")
score2 <- read.csv("# diseases_pca_score.csv")
# Adding indegree to score dataset
dat2 <- merge(score2, Centrality, by.x= "FHA", by.y="ID")

# Subsetting only to SW sites
d2 <- dat2[dat2$Type=="SW",]
# Remove variables with missing values (otherwiese stan models crash) and variables not used
d2 <- d2[,c("score2","No_inf", "No.outbreaks", "indeg", "windeg")]

# Centering
# construct centered predictor
d2$score2_c <- d2$score2 - mean(d2$score2)
d2$windeg_c <- d2$windeg - mean(d2$windeg)
d2$windeg_cat <- ifelse(d$windeg>1,1,0)

# Changing name of No.outbreaks variable, as the dot annoys stan
colnames(d2)[3] <- "No_outbreaks"

## re-estimate
# Full model
m8.1 <- map2stan(
  alist(No_outbreaks ~ dpois( lambda ),
        log(lambda) <- a + bs*score2_c+ bi*windeg_cat + bsi*score2_c*windeg_cat,
        a ~ dnorm(0,100),
        c(bs,bi,bsi) ~ dnorm(0,1)
  ),
  data=d2, iter=3000 , warmup=1000 , chains=4 ,
  control = list(adapt_delta = 0.80))

precis(m8.1, prob=0.95, digits=3)
pairs(m5.1)

# No interaction
m8.2 <- map2stan(
  alist(No_outbreaks ~ dpois( lambda ),
        log(lambda) <- a + bs*score2_c+ bi*windeg_cat ,
        a ~ dnorm(0,100),
        c(bs,bi) ~ dnorm(0,1)
  ),
  data=d2, iter=3000 , warmup=1000 , chains=4 ,
  control = list(adapt_delta = 0.80))

precis(m8.2, prob=0.95, digits=3)
pairs(m8.2)


# No windeg
m8.3 <- map2stan(
  alist(No_outbreaks ~ dpois( lambda ),
        log(lambda) <- a + bs*score2_c,
        a ~ dnorm(0,100),
        c(bs) ~ dnorm(0,1)
  ),
  data=d2, iter=3000 , warmup=1000 , chains=4 )

precis(m8.3, prob=0.95, digits=3)
pairs(m8.3)

# No biosec
m8.4 <- map2stan(
  alist(No_outbreaks ~ dpois( lambda ),
        log(lambda) <- a + bi*windeg_cat ,
        a ~ dnorm(0,100),
        c(bi) ~ dnorm(0,1)
  ),
  data=d2, iter=3000 , warmup=1000 , chains=4 ,
  control = list(adapt_delta = 0.80))

precis(m8.4, prob=0.95, digits=3)
pairs(m8.4)

# Null model
m5.5 <- map2stan(
  alist(No_outbreaks ~ dpois( lambda ),
        log(lambda) <- a ,
        a ~ dnorm(0,100)
  ),
  data=d2, iter=3000 , warmup=1000 , chains=4 ,
  control = list(adapt_delta = 0.80))

precis(m5.5, prob=0.95, digits=3)
plot(m5.5)


# Comparing models

diseases.compare.stan_pca <- compare(m8.1,m8.2,m8.3,m8.4,m5.5,n=1e4)
plot(diseases.compare.stan_pca)

# bis a bis models with score based on full data and based on partial data
diseases.compare.stan_pca_vs_full <- compare(m8.1,m8.2,m8.3,m8.4,m6.1, m6.2, m5.3, m6.4, m5.5,n=1e4)
plot(diseases.compare.stan_pca_vs_full)


## Comparing high and low biosec for high contact site
post8.2 <- extract.samples(m8.2)
dis_high_bs <- exp(post8.2$a + post8.2$bi+ post8.2$bs*4.4560)
dis_low_bs <- exp(post8.2$a + post8.2$bi+ post8.2$bs*-1.0320)

# Difference in number of diseases
dis_diff <- dis_low_bs - dis_high_bs
sum(dis_diff > 0)/length(dis_diff)
plot(density(dis_diff), xlab="Difference", main="")
abline(v=0, lty=2)

## Comparing high and low biosec for low contact site
dis_high_bs_lowc <- exp(post8.2$a + post8.2$bs*4.4560)
dis_low_bs_lowc <- exp(post8.2$a + post8.2$bs*-1.0320)
# Difference in number of diseases
dis_diff_lowc <- dis_low_bs_lowc - dis_high_bs_lowc
sum(dis_diff_lowc > 0)/length(dis_diff_lowc)
plot(density(dis_diff_lowc), xlab="Difference", main="")
abline(v=0, lty=2)


## Comparing high and low contact for high biosec
dis_high_bi <- exp(post8.2$a + post8.2$bi+ post8.2$bs*4.4560)
dis_low_bi <- exp(post8.2$a + post8.2$bs*4.4560)
# Difference in number of diseases
dis_diff_bi <- dis_high_bi - dis_low_bi
sum(dis_diff_bi > 0)/length(dis_diff_bi)
plot(density(dis_diff_bi), xlab="Difference", main="")
abline(v=0, lty=2)

## Comparing high and low contact for low biosec
dis_high_bi_lowbs <- exp(post8.2$a + post8.2$bi+ post8.2$bs*-1.0320)
dis_low_bi_lowbs <- exp(post8.2$a + post8.2$bs*-1.0320)
# Difference in number of diseases
dis_diff_bi_lowbs <- dis_high_bi_lowbs - dis_low_bi_lowbs
sum(dis_diff_bi_lowbs > 0)/length(dis_diff_bi_lowbs)
plot(density(dis_diff_bi_lowbs), xlab="Difference", main="")
abline(v=0, lty=2)

# Plotting model predictions
pch <- ifelse( d$windeg > 1 , 16 , 1 )
plot( d2$score2_c , d2$No_inf , col=rangi2 , pch=pch ,
      xlab="Biosecurity score" , ylab="No of different diseases", ylim=c(0, 11) )

# sequence of biosecurity scores and contact rate to compute over
biosec.seq <- seq( from=-11.5 , to=5.5 , length.out=80 )

# compute trend
# High contact rate
d.pred <- data.frame(
  score2_c = biosec.seq,
  windeg_cat = 1
)

lambda.pred.h <- link( m8.2, data=d.pred)
lambda.med <- apply( lambda.pred.h , 2 , median )
lambda.PI <- apply( lambda.pred.h , 2 , HPDI )

# plot predicted trend
lines( biosec.seq , lambda.med , col=rangi2 )
shade( lambda.PI , biosec.seq , col=col.alpha(rangi2,0.2) )

# Low contact rate
d.pred.cat_l <- data.frame(
  score2_c = biosec.seq,
  windeg_cat = 0
)

lambda.pred.l <- link(m8.2, data=d.pred.cat_l)
lambda.med.l <- apply( lambda.pred.l, 2 , median )
lambda.PI.l <- apply( lambda.pred.l, 2 , HPDI )

lines( biosec.seq , lambda.med.l , lty=2 )
shade( lambda.PI.l , biosec.seq , col=col.alpha("black",0.1) )



# Multi-level adjusting for manager ID
# Full model...similar results, but poorer fitting
# Perhaps due to over fitting: 15 parameters with just 17 obs
d$Manager <- as.integer(d$Manager)

m8.1 <- map2stan(
  alist(No_outbreaks ~ dpois( lambda ),
        log(lambda) <- a + a_manager[Manager]+ bs*Score_surv100_c+ bi*windeg_c + bsi*Score_surv100_c*windeg_c,
        a_manager[Manager] ~dnorm(0, sigma),
        a ~ dnorm(0,100),
        c(bs,bi,bsi) ~ dnorm(0,1),
        sigma ~ dcauchy(0, 1)
  ),
  data=d, iter=20000 , warmup=1000 , chains=1 ,start=list(a=0, bs=0, bi=0, bsi=0),
  control = list(adapt_delta = 0.90))

precis(m8.1, prob=0.95, digits=3, depth=2)
pairs(m8.1)

# No interaction
m8.2 <- map2stan(
  alist(No_outbreaks ~ dpois( lambda ),
        log(lambda) <- a + a_manager[Manager]+ bs*Score_surv100_c+ bi*windeg_c ,
        a_manager[Manager] ~dnorm(0, sigma),
        a ~ dnorm(0,100),
        c(bs,bi) ~ dnorm(0,1),
        sigma ~ dcauchy(0, 2.5)
  ),
  data=d, iter=20000 , warmup=1000 , chains=1 ,start=list(a=0, bs=0, bi=0),
  control = list(adapt_delta = 0.90))

precis(m8.2, prob=0.95, digits=3, depth=2)
pairs(m8.2)

# No indeg
m8.3 <- map2stan(
  alist(No_outbreaks ~ dpois( lambda ),
        log(lambda) <- a + a_manager[Manager]+ bs*Score_surv100_c ,
        a_manager[Manager] ~dnorm(0, sigma),
        a ~ dnorm(0,100),
        c(bs) ~ dnorm(0,1),
        sigma ~ dcauchy(0, 2.5)
  ),
  data=d, iter=20000 , warmup=1000 , chains=1 ,start=list(a=0, bs=0),
  control = list(adapt_delta = 0.90))

precis(m8.3, prob=0.95, digits=3, depth=2)
pairs(m8.3)

# No biosec
m8.4 <- map2stan(
  alist(No_outbreaks ~ dpois( lambda ),
        log(lambda) <- a + a_manager[Manager]+ bi*windeg_c ,
        a_manager[Manager] ~dnorm(0, sigma),
        a ~ dnorm(0,100),
        c(bi) ~ dnorm(0,1),
        sigma ~ dcauchy(0, 2.5)
  ),
  data=d, iter=20000 , warmup=1000 , chains=1 ,start=list(a=0, bi=0),
  control = list(adapt_delta = 0.90))

precis(m8.4, prob=0.95, digits=3, depth=2)
pairs(m8.4)


# Null model
m8.5 <- map2stan(
  alist(No_outbreaks ~ dpois( lambda ),
        log(lambda) <- a + a_manager[Manager],
        a_manager[Manager] ~dnorm(0, sigma),
        a ~ dnorm(0,100),
        sigma ~ dcauchy(0, 2.5)
  ),
  data=d, iter=20000 , warmup=1000 , chains=1 ,start=list(a=0),
  control = list(adapt_delta = 0.90))

precis(m8.5, prob=0.95, digits=3, depth=2)
pairs(m8.5)

# Comparing with non multi-level ones. The latter are always better
comp <- compare(m8.1, m8.2, m8.3, m8.4, m8.5,m5.1, m5.2, m5.3, m5.4, m5.5)
compare(m5.1, m8.1)
compare(m5.2, m8.2)
compare(m5.3, m8.3)
compare(m5.4, m8.4)
compare(m5.5, m8.5)


# Running a generalized poisson model in Rstan
### Ben's advice

# Data
intercept <- rep(1, nrow(d))
inter <- d$Score_surv100_c*d$windeg_c
covs <- as.matrix(data.frame(intercept, as.vector(d$Score_surv100_c), d$windeg_c,inter) )
X <- covs                  
y <- d$No_outbreaks 

dat <- list(X=X, y=y, N=nrow(X), K=ncol(X))


write.table(dat$X, "covariates.txt", row.names=F)
write.table(dat$y, "outcome.txt", row.names=T)
write.table(dat$N, "N.txt", row.names=F)
write.table(dat$K, "K.txt", row.names=F)


# Use this data in the meantime, as I don;t want to re-run the whole code
intercept <- rep(1, 17)
Score_surv100_c=c(-2.7588458, -2.7588458,  3.4339806,  3.0432020,  2.7448574,  3.7356885, -1.7610177, -4.1047376, -2.4006051, -2.4006051, -5.1737057, -0.7110643, -0.3697673, -0.7110643,  3.5651855,  3.9328326, 2.6945123)
windeg_c=c(1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1)
inter <- Score_surv100_c*windeg_c
covs <- as.matrix(data.frame(intercept, Score_surv100_c, windeg_c, inter) )
X <- covs                  
y <- c(1, 1, 0, 2, 1, 0, 2, 2, 1, 2, 3, 1, 1, 1, 2, 2, 2) 

dat <- list(X=X, y=as.integer(y), N=17, K=4, m=5)


#######################################################################################
################# It should be fitted as a multi-level model, with farm random effects, correlated with distance!
library(shinystan)
library(rstan)


data <- list(score=as.vector(d$Score_surv100_c), windeg=d$windeg_c, inter=as.vector(d$inter), N=nrow(d), 
             y=d$No_outbreaks, site = seq(1:17), n_site=nrow(d))
data.pca <- list(score=as.vector(d$Score_55_c), windeg=d$windeg_c, inter=as.vector(d$Score_55_c*d$windeg_c), N=nrow(d), 
             y=d$No_outbreaks, site = seq(1:17), n_site=nrow(d))


data2 <- list(score=as.vector(d2$Score_surv100_c), windeg=d2$windeg_c, inter=as.vector(d2$inter), N=nrow(d2), 
              y=d2$No_outbreaks, site = seq(1:17), n_site=nrow(d2))

data3 <- list(score=as.vector(d3$Score_surv100_c), windeg=d3$windeg_c, inter=as.vector(d3$inter), N=nrow(d3), 
              y=d3$No_outbreaks, site = seq(1:34), n_site=nrow(d3), type=as.integer(d3$Type),
              N_type=length(unique(d3$Type)))


stancode<- "

data {
int<lower=0> N;
vector[N] score;
vector[N] windeg;
vector[N] inter;
vector<lower=0>[N] y;
int<lower=1> site[N];
int <lower=0> N_type;
int type[N];
}
transformed data {
vector[N] ym1 = y - 1;
}
parameters {
real a;
real bs;
real bw;
real bsw;
vector[N_type] a_site;
vector[N_type] bs_site;
vector[N_type] bw_site;
vector[N_type] bsw_site;
vector<lower=0>[4] sigma_site;
corr_matrix[4] Rho;
real<lower=-1, upper=1> lambda[N_type];
}
transformed parameters{
vector[4] v_a_sitecoef_site[N_type];
vector[4] Mu_abm;
cov_matrix[4] SRS_sigma_siteRho;
    for ( j in 1:N_type ) {
v_a_sitecoef_site[j,1] = a_site[j];
v_a_sitecoef_site[j,2] = bs_site[j];
v_a_sitecoef_site[j,3] = bw_site[j];
v_a_sitecoef_site[j,4] = bsw_site[j];
    }

for ( j in 1:4 ) {
        Mu_abm[1] = a;
        Mu_abm[2] = bs;
        Mu_abm[3] = bw;
        Mu_abm[4] = bsw;
    }

SRS_sigma_siteRho = quad_form_diag(Rho,sigma_site);
}


model {
vector[N] eta;
vector[N] theta;
vector[N] theta_lambda_y;
vector[N_type] A;
vector[N_type] BS;
vector[N_type] BW;
vector[N_type] BSW;
vector[N] log_1m_lambda;

for (i in 1:N){
A[type[i]] = a + a_site[type[i]];
BS[type[i]] = bs + bs_site[type[i]];
BW[type[i]] = bw + bw_site[type[i]];
BSW[type[i]] = bsw + bsw_site[type[i]];
eta[i]= A[type[i]] + BS[type[i]]*score[i] + BW[type[i]]*windeg[i] + BSW[type[i]]*inter[i];
theta[i] = exp(eta[i]) * (1 - lambda[type[i]]);
theta_lambda_y[i] = theta[i] + lambda[type[i]] * y[i];
log_1m_lambda[i] = log1m(lambda[type[i]]);
}


target += (eta + log_1m_lambda + ym1 .* log(theta_lambda_y) - theta_lambda_y);

lambda ~ normal(0,1);
a ~ normal(0,0.5);
sigma_site ~ cauchy(0,1);
bs ~ normal(0,0.5);
bw ~ normal(0,0.5);
bsw ~ normal(0,0.5);
Rho ~ lkj_corr( 3 );
v_a_sitecoef_site ~multi_normal(Mu_abm,SRS_sigma_siteRho);

}

generated quantities{
vector[N_type] A;
vector[N_type] BS;
vector[N_type] BW;
vector[N_type] BSW;
vector[N] log_lik;
vector[N] theta;
vector[N] eta;

for (i in 1:N){
A[type[i]] = a + a_site[type[i]];
BS[type[i]] = bs + bs_site[type[i]];
BW[type[i]] = bw + bw_site[type[i]];
BSW[type[i]] = bsw + bsw_site[type[i]];
eta[i]= A[type[i]] + BS[type[i]]*score[i] + BW[type[i]]*windeg[i] + BSW[type[i]]*inter[i];
log_lik[i] = (eta[i] + log1m(lambda[type[i]]) + ym1[i] .* log((exp(eta[i]) * (1 - lambda[type[i]])) + lambda[type[i]] * y[i]) - ((exp(eta[i]) * (1 - lambda[type[i]])) + lambda[type[i]] * y[i]) - lgamma(y[i] + 1));
theta[i] = exp(eta[i]) * (1 - lambda[type[i]]);
  }
}

"

stancode2 <- "
functions{
real lb(vector score, vector windeg, vector inter,  real intercept, real Score, real NoProviders, real Interaction){
real max_mu = max(exp(intercept + Score*score +NoProviders*windeg + Interaction*inter));
if (max_mu <= 1.5) return max_mu / (max_mu - 3);
return -1;    
}
}
  data {
int<lower=0> N;
vector[N] score;
vector[N] windeg;
vector[N] inter;
vector<lower=0>[N] y;
}
transformed data {
vector[N] ym1 = y - 1;
}
parameters {
real intercept;
real Score;
real NoProviders;
real Interaction;
real <lower=lb(score, windeg, inter, intercept, Score, NoProviders, Interaction), upper = 1> lambda;
}
model {
vector[N] eta;
vector[N] theta;
vector[N] theta_lambda_y;


for (i in 1:N){
eta[i]= intercept + Score*score[i] +NoProviders*windeg[i] + Interaction*inter[i];
}
theta = exp(eta) * (1 - lambda);
theta_lambda_y = theta + lambda * y;

target += (eta + log1m(lambda) + ym1 .* log(theta_lambda_y) - theta_lambda_y);

lambda ~ normal (0,1);
intercept ~ normal(0,0.5);
Score ~ normal(0,0.5);
NoProviders ~ normal(0,0.5);
Interaction ~ normal(0,0.5);
}

generated quantities{
vector[N] log_lik;
vector[N] theta;
vector[N] eta;

for (i in 1:N){

eta[i]= intercept + Score*score[i] +NoProviders*windeg[i] + Interaction*inter[i];
log_lik[i] = (eta[i] + log1m(lambda) + ym1[i] .* log((exp(eta[i]) * (1 - lambda)) + lambda * y[i]) - ((exp(eta[i]) * (1 - lambda)) + lambda * y[i]) - lgamma(y[i] + 1));
theta[i] = exp(eta[i]) * (1 - lambda);
}
}
"


stancode3 <- "
functions{
real lb(vector score, vector windeg, real a, real bs, real bw){
real max_mu = max(exp(a + bs*score +bw*windeg));
if (max_mu <= 1.5) return max_mu / (max_mu - 3);
return -1;    
}
}  
  data {
  int<lower=0> N;
  vector[N] score;
  vector[N] windeg;
  vector[N] inter;
  vector<lower=0>[N] y;
}
transformed data {
  vector[N] ym1 = y - 1;
}
parameters {
  real a;
  real bs;
  real bw;
  real <lower=lb(score, windeg, a, bs, bw), upper = 1> lambda;
}
model {
  vector[N] eta;
  vector[N] theta;
  vector[N] theta_lambda_y;
 
  
  for (i in 1:N){
    eta[i]= a + bs*score[i] +bw*windeg[i];
  }
  theta = exp(eta) * (1 - lambda);
  theta_lambda_y = theta + lambda * y;
  
  target += (eta + log1m(lambda) + ym1 .* log(theta_lambda_y) - theta_lambda_y);
  
  lambda ~ normal(0,1);
a ~ normal(0,0.5);
bs ~ normal(0,0.5);
bw ~ normal(0,0.5);
}

generated quantities{
  vector[N] log_lik;
  vector[N] theta;
  vector[N] eta;
  
  for (i in 1:N){
    
    eta[i]= a + bs*score[i] +bw*windeg[i];
    log_lik[i] = (eta[i] + log1m(lambda) + ym1[i] .* log((exp(eta[i]) * (1 - lambda)) + lambda * y[i]) - ((exp(eta[i]) * (1 - lambda)) + lambda * y[i]) - lgamma(y[i] + 1));
    theta[i] = exp(eta[i]) * (1 - lambda);
  }
}
"

stancode4 <- "
functions{
real lb(vector score, real a, real bs){
real max_mu = max(exp(a + bs*score));
if (max_mu <= 1.5) return max_mu / (max_mu - 3);
return -1;    
}
}  
  data {
  int<lower=0> N;
  vector[N] score;
  vector[N] inter;
  vector<lower=0>[N] y;
}
transformed data {
  vector[N] ym1 = y - 1;
}
parameters {
  real a;
  real bs;
  real<lower=lb(score, a, bs), upper = 1> lambda;
}
model {
  vector[N] eta;
  vector[N] theta;
  vector[N] theta_lambda_y;
 
  
  for (i in 1:N){
    eta[i]= a + bs*score[i];
  }
  theta = exp(eta) * (1 - lambda);
  theta_lambda_y = theta + lambda * y;
  
  target += (eta + log1m(lambda) + ym1 .* log(theta_lambda_y) - theta_lambda_y);
  
  lambda ~ normal(0,1);
a ~ normal(0,0.5);
bs ~ normal(0,0.5);
}

generated quantities{
  vector[N] log_lik;
  vector[N] theta;
  vector[N] eta;
  
  for (i in 1:N){
    
    eta[i]= a + bs*score[i];
    log_lik[i] = (eta[i] + log1m(lambda) + ym1[i] .* log((exp(eta[i]) * (1 - lambda)) + lambda * y[i]) - ((exp(eta[i]) * (1 - lambda)) + lambda * y[i]) - lgamma(y[i] + 1));
    theta[i] = exp(eta[i]) * (1 - lambda);
  }
}
"

stancode5 <- "
functions{
real lb(vector windeg, real a, real bw){
real max_mu = max(exp(a + bw*windeg));
if (max_mu <= 1.5) return max_mu / (max_mu - 3);
return -1;    
}
}  
  data {
  int<lower=0> N;
  vector[N] windeg;
  vector<lower=0>[N] y;
}
transformed data {
  vector[N] ym1 = y - 1;
}
parameters {
  real a;
  real bw;
  real <lower=lb(windeg, a, bw), upper = 1> lambda;
}
model {
  vector[N] eta;
  vector[N] theta;
  vector[N] theta_lambda_y;
 
  
  for (i in 1:N){
    eta[i]= a + bw*windeg[i];
  }
  theta = exp(eta) * (1 - lambda);
  theta_lambda_y = theta + lambda * y;
  
  target += (eta + log1m(lambda) + ym1 .* log(theta_lambda_y) - theta_lambda_y);
  
  lambda ~ normal(0,1);
a ~ normal(0,0.5);
bw ~ normal(0,0.5);
}

generated quantities{
  vector[N] log_lik;
  vector[N] theta;
  vector[N] eta;
  
  for (i in 1:N){
    
    eta[i]= a + bw*windeg[i];
    log_lik[i] = (eta[i] + log1m(lambda) + ym1[i] .* log((exp(eta[i]) * (1 - lambda)) + lambda * y[i]) - ((exp(eta[i]) * (1 - lambda)) + lambda * y[i]) - lgamma(y[i] + 1));
    theta[i] = exp(eta[i]) * (1 - lambda);
  }
}
"

stancode6 <- "
functions{
real lb(real a){
real max_mu = exp(a);
if (max_mu <= 1.5) return max_mu / (max_mu - 3);
return -1;    
}
}
data {
  int<lower=0> N;
  vector<lower=0>[N] y;
}
transformed data {
  vector[N] ym1 = y - 1;
}
parameters {
  real a;
  real <lower=lb(a), upper = 1> lambda;
}
model {
  vector[N] eta;
  vector[N] theta;
  vector[N] theta_lambda_y;
 
  
  for (i in 1:N){
    eta[i]= a ;
  }
  theta = exp(eta) * (1 - lambda);
  theta_lambda_y = theta + lambda * y;
  
  target += (eta + log1m(lambda) + ym1 .* log(theta_lambda_y) - theta_lambda_y);
  
lambda ~ normal(0,1);
a ~ normal(0,0.5);
}

generated quantities{
  vector[N] log_lik;
  vector[N] theta;
  vector[N] eta;
  
  for (i in 1:N){
    
    eta[i]= a ;
    log_lik[i] = (eta[i] + log1m(lambda) + ym1[i] .* log((exp(eta[i]) * (1 - lambda)) + lambda * y[i]) - ((exp(eta[i]) * (1 - lambda)) + lambda * y[i]) - lgamma(y[i] + 1));
    theta[i] = exp(eta[i]) * (1 - lambda);
  }
}
"
stancode7 <- 
"data {
  int<lower=0> N;
vector[N] score;
vector[N] windeg;
vector[N] inter;
int y [N];
}
parameters{
    real a;
    real bs;
    real bw;
    real bsw;
}
model{
    vector[N] lambda;
    bsw ~ normal( 0 , 1 );
    bw ~ normal( 0 , 1 );
    bs ~ normal( 0 , 1 );
    a ~ normal( 0 , 1 );
    for ( i in 1:N ) {
        lambda[i] = a + bs * score[i] + bw * windeg[i] + bsw *inter[i];
    }
    y ~ poisson_log( lambda );
}
generated quantities{
    vector[N] lambda;
    vector[N] log_lik;
    for ( i in 1:N ) {
        lambda[i] = a + bs * score[i] + bw * windeg[i] + bsw *inter[i];
        log_lik[i] = poisson_log_lpmf( y[i] | lambda[i] );
    }
}"


n_chains <- 4
start <- list(lambda=0.0, intercept=0, Score=0, NoProviders=0, Interaction=0)
init <- list()
for ( i in 1:n_chains ) init[[i]] <- start

start <- list(sigma_site=rep(1,4), lambda=0.0, a=0, bs=0, bw=0, bsw=0)

m9.1 <- stan(model_code = stancode, data=data3, iter=4000, chains=4,cores=4, thin=1,
             warmup=2000, control = list(adapt_delta = 0.9, stepsize=0.5))


m9.2 <- stan(model_code = stancode2, data=data, iter=16000, chains=1,cores=1, thin=1,
             warmup=8000, control = list(adapt_delta = 0.9, stepsize=0.5), init=init)

Full_model <- m9.2
Main_effects <- stan(model_code = stancode3, data=data, iter=4000, chains=4,cores=4, thin=1,
             warmup=2000, control = list(adapt_delta = 0.9, stepsize=0.5), init=init)
Score <- stan(model_code = stancode4, data=data, iter=4000, chains=4,cores=4, thin=1,
             warmup=2000, control = list(adapt_delta = 0.99, stepsize=0.5), init=init)
Windeg <- stan(model_code = stancode5, data=data, iter=4000, chains=4,cores=4, thin=1,
             warmup=2000, control = list(adapt_delta = 0.9, stepsize=0.5), init=init)

Null_model <- stan(model_code = stancode6, data=data, iter=4000, chains=4,cores=4, thin=1,
               warmup=2000, control = list(adapt_delta = 0.99, stepsize=0.5), init=init)

mpoiss <- stan(model_code = stancode7, data=data, iter=4000, chains=4,cores=4, thin=1,
               warmup=2000, control = list(adapt_delta = 0.9, stepsize=1))

m9.2.pca <- stan(model_code = stancode2, data=data.pca, iter=4000, chains=4,cores=4, thin=1,
             warmup=2000, control = list(adapt_delta = 0.9, stepsize=0.5), init=init)


# Models for not seawater sites

m9.2.fw <- stan(model_code = stancode2, data=data2, iter=16000, chains=1,cores=1, thin=1,
             warmup=8000, control = list(adapt_delta = 0.9, stepsize=0.5), init=init)

Full_model.fw <- m9.2
Main_effects.fw <- stan(model_code = stancode3, data=data2, iter=4000, chains=4,cores=4, thin=1,
                     warmup=2000, control = list(adapt_delta = 0.9, stepsize=0.5), init=init)
Score.fw <- stan(model_code = stancode4, data=data2, iter=4000, chains=4,cores=4, thin=1,
              warmup=2000, control = list(adapt_delta = 0.9, stepsize=0.5), init=init)
Windeg.fw <- stan(model_code = stancode5, data=data2, iter=4000, chains=4,cores=4, thin=1,
               warmup=2000, control = list(adapt_delta = 0.9, stepsize=0.5), init=init)

Null_model.fw <- stan(model_code = stancode6, data=data2, iter=4000, chains=4,cores=4, thin=1,
                   warmup=2000, control = list(adapt_delta = 0.9, stepsize=0.5), init=init)

mpoiss.fw <- stan(model_code = stancode7, data=data2, iter=4000, chains=4,cores=4, thin=1,
               warmup=2000, control = list(adapt_delta = 0.9, stepsize=1))

print(m9.1, probs=c(0.025, 0.975), pars=c("a", "bs_site", "bw_site", "bsw_site", "A", "BS", "BW", "BSW","lambda"))
print(m9.2, probs=c(0.025, 0.975), pars=c("intercept", "Score", "NoProviders", "Interaction", "lambda"))
print(Main_effects, probs=c(0.025, 0.975), pars=c("log_lik", "theta", "eta", "log-posterior"), include=F)
print(Score, probs=c(0.025, 0.975), pars=c("log_lik", "theta", "eta", "log-posterior"), include=F)
print(Windeg, probs=c(0.025, 0.975), pars=c("log_lik", "theta", "eta", "log-posterior"), include=F)
print(Null_model, probs=c(0.025, 0.975), pars=c("log_lik", "theta", "eta", "log-posterior"), include=F)
print(m9.2.pca, probs=c(0.025, 0.975), pars=c("intercept", "Score", "NoProviders", "Interaction", "lambda"))


pairs(m9.2.fw, pars=c("intercept", "Score", "NoProviders", "Interaction", "lambda"))

traceplot(m9.1, inc_warmup=F, pars=c("log_lik", "theta", "eta", "log-posterior","sigma_site", "Rho"), include=F)
traceplot(m9.2, pars=c("log_lik", "theta", "eta", "log-posterior"), include=F)
traceplot(Main_effects, pars=c("log_lik", "theta", "eta", "log-posterior"), include=F)
traceplot(Score, pars=c("log_lik", "theta", "eta", "log-posterior"), include=F)
traceplot(Windeg, pars=c("log_lik", "theta", "eta", "log-posterior"), include=F)
traceplot(Null_model, pars=c("log_lik", "theta", "eta", "log-posterior"), include=F)
traceplot(m9.2.pca, pars=c("log_lik", "theta", "eta", "log-posterior"), include=F)

print(m9.2.fw, probs=c(0.025, 0.975), pars=c("intercept", "Score", "NoProviders", "Interaction", "lambda"))
print(Main_effects.fw, probs=c(0.025, 0.975), pars=c("log_lik", "theta", "eta", "log-posterior"), include=F)
print(Score.fw, probs=c(0.025, 0.975), pars=c("log_lik", "theta", "eta", "log-posterior"), include=F)
print(Windeg.fw, probs=c(0.025, 0.975), pars=c("log_lik", "theta", "eta", "log-posterior"), include=F)
print(Null_model.fw, probs=c(0.025, 0.975), pars=c("log_lik", "theta", "eta", "log-posterior"), include=F)

traceplot(m9.2.fw, pars=c("log_lik", "theta", "eta", "log-posterior"), include=F)
traceplot(Main_effects.fw, pars=c("log_lik", "theta", "eta", "log-posterior"), include=F)
traceplot(Score.fw, pars=c("log_lik", "theta", "eta", "log-posterior"), include=F)
traceplot(Windeg.fw, pars=c("log_lik", "theta", "eta", "log-posterior"), include=F)
traceplot(Null_model.fw, pars=c("log_lik", "theta", "eta", "log-posterior"), include=F)


stan_plot(m9.2, pars=c("intercept", "Score", "NoProviders", "Interaction", "lambda"))
stan_plot(m9.2.fw, pars=c("intercept", "Score", "NoProviders", "Interaction", "lambda"))
stan_plot(m9.1, pars=c("bs_site", "bw_site", "bsw_site", "lambda"))
stan_plot(m9.2.pca, pars=c("intercept", "Score", "NoProviders", "Interaction", "lambda"))
stan_par(m9.2, par="bs")
stan_diag(m9.2)
stan_rhat(m9.2)
stan_ess(m9.2)
stan_mcse(m9.2)



# Estimating WAIC from posterior samples
library(loo)
log_lik1 <- extract_log_lik(m9.2)
(Full_loo <- loo(log_lik1))
log_lik2 <- extract_log_lik(Main_effects)
(Main_loo <- loo(log_lik2))
log_lik3 <- extract_log_lik(Score)
(Score_loo <- loo(log_lik3))
log_lik4 <- extract_log_lik(Windeg)
(Windeg_loo <- loo(log_lik4))
log_lik5 <- extract_log_lik(Null_model)
(Null_loo <- loo(log_lik5))
log_lik6 <- extract_log_lik(mpoiss)
(mpoiss_loo <- loo(log_lik6))

compare(Full_loo, Main_loo, Score_loo, Windeg_loo, Null_loo)


log_lik1 <- extract_log_lik(m9.2.fw)
(Full_loo <- loo(log_lik1))
log_lik2 <- extract_log_lik(Main_effects.fw)
(Main_loo <- loo(log_lik2))
log_lik3 <- extract_log_lik(Score.fw)
(Score_loo <- loo(log_lik3))
log_lik4 <- extract_log_lik(Windeg.fw)
(Windeg_loo <- loo(log_lik4))
log_lik5 <- extract_log_lik(Null_model.fw)
(Null_loo <- loo(log_lik5))
log_lik6 <- extract_log_lik(mpoiss.fw)
(mpoiss_loo <- loo(log_lik6))

compare(Full_loo, Main_loo, Score_loo, Windeg_loo, Null_loo)


shiny_m9.2 <- as.shinystan(m9.2)
shiny_pois <- as.shinystan(mpoiss)
launch_shinystan(shiny_m9.2)

post <- extract(m9.2)
post.fw <- extract(m9.2.fw)

# MAKE A DOT PLOT OF POSTERIOR DISTRIBUTION OF PARAMETER ESTIMATES
# SW model
param <- as.data.frame(post); param <- param[,1:5]
colnames(param) <- c("Intercept", "Score", "No Suppliers", "Interaction", "Lambda")
var.labels <- names(param)
med <- apply(param, 2, median)
pi.95 <- apply(param, 2, PI, prob=0.95)
lb <- pi.95[1,]; ub <- pi.95[2,]
pi.80 <- apply(param, 2, PI, prob=0.66)
lb.80 <- pi.80[1,]; ub.80 <- pi.80[2,]

sum.param <- data.frame(factor(var.labels, levels=names(param)),med, lb, ub, lb.80, ub.80)
colnames(sum.param) <- c("param", "Estimate", "lb", "ub", "lb.80", "ub.80")
row.names(sum.param) <- NULL


# FW model
param.fw <- as.data.frame(post.fw); param.fw <- param.fw[,1:5]
colnames(param.fw) <- c("Intercept", "Score", "No Suppliers", "Interaction", "Lambda")
var.labels.fw <- names(param.fw)
med.fw <- apply(param.fw, 2, median)
pi.95.fw <- apply(param.fw, 2, PI, prob=0.95)
lb.fw <- pi.95.fw[1,]; ub.fw <- pi.95.fw[2,]
pi.80.fw <- apply(param.fw, 2, PI, prob=0.66)
lb.80.fw <- pi.80.fw[1,]; ub.80.fw <- pi.80.fw[2,]

sum.param.fw <- data.frame(factor(var.labels.fw, levels=names(param.fw)),med.fw, lb.fw, ub.fw, lb.80.fw, ub.80.fw)
colnames(sum.param.fw) <- c("param", "Estimate", "lb", "ub", "lb.80", "ub.80")
row.names(sum.param.fw) <- NULL

# Plot

plot1 <- dotplot(param~Estimate, data=sum.param, xlim=c(-max(sum.param$ub)-0.02, 
        max(sum.param$ub)+0.02), scales=list(y=list(cex=8/12), x=list(cex=8/12)), xlab="", 
        panel=function(x,y){
  panel.xyplot(x, y, pch=16, cex=1, col="black", grid=T)
  panel.segments(sum.param$lb, as.numeric(y), sum.param$ub, as.numeric(y), lty=1, 
                 col="black", lwd=1)
  panel.segments(sum.param$lb.80, as.numeric(y), sum.param$ub.80, as.numeric(y), lty=1, 
                 lwd=2, col="red3")
  panel.xyplot(x, y, pch=16, cex=1, col="black")
  panel.abline(v=0, col=1, lty=2)
})

plot2 <- dotplot(param~Estimate, data=sum.param.fw, xlim=c(-max(sum.param.fw$ub)-0.02, 
        max(sum.param.fw$ub)+0.02), scales=list(y=list(cex=8/12), x=list(cex=8/12)), xlab="", 
        panel=function(x,y){
  panel.xyplot(x, y, pch=16, cex=1, col="black", grid=T)
  panel.segments(sum.param.fw$lb, as.numeric(y), sum.param.fw$ub, as.numeric(y), lty=1, 
                 col="black", , lwd=1)
  panel.segments(sum.param.fw$lb.80, as.numeric(y), sum.param.fw$ub.80, as.numeric(y), lty=1, 
                 lwd=2, col="red3")
  panel.xyplot(x, y, pch=16, cex=1, col="black")
  panel.abline(v=0, col=1, lty=2)
})

tiff(filename = "Fig1.tiff",
               width = 7.5, height = 3, units = "in", pointsize = 12,
               compression = "lzw",
               bg = "white", res = 300, type = "windows")

print(plot1, split = c(1,1,2,1), more=T)
print(plot2, split = c(2,1,2,1), more=F)

dev.off()

# Plotting prior and posterior distribution
Intercept <- post$intercept
Score <- post$Score
NoProviders <- post$NoProviders
Interaction <- post$Interaction
lambda <- post$lambda
theta <- post$theta

Intercept.fw <- post.fw$intercept
Score.fw <- post.fw$Score
NoProviders.fw <- post.fw$NoProviders
Interaction.fw <- post.fw$Interaction
lambda.fw <- post.fw$lambda
theta.fw <- post.fw$theta

par(mfrow=c(2,2))
prior_eta <- rnorm(n=10000, mean=0, sd=0.5)
prior_lambda <- rnorm(n=10000, mean=0, sd=1)

plot(density(prior_eta), ylim=c(0,3), main="Intercept")
lines(density(Intercept), col="blue")
lines(density(Intercept.fw), col="red")

plot(density(prior_eta), ylim=c(0,2), main="Score")
lines(density(Score), col="blue")
lines(density(Score.fw), col="red")


plot(density(prior_eta), ylim=c(0,2), main="NoProviders")
lines(density(NoProviders), col="blue")
lines(density(NoProviders.fw), col="red")

plot(density(prior_eta), ylim=c(0,2), main="Interaction")
lines(density(Interaction), col="blue")
lines(density(Interaction.fw), col="red")

plot(density(prior_lambda), ylim=c(0,2), main="lambda")
lines(density(lambda), col="blue")
lines(density(lambda.fw), col="red")

## Simulate counts
library(RMKdiscrete)

kk <- matrix(rep(NA,8000*17), ncol=17)

for (i in 1:17){   
    kk[,i] <- rLGP(1,theta[,i],lambda)
  }

tiff(filename = "Fig4.tiff",
               width = 7.5, height = 4, units = "in", pointsize = 12,
               compression = "lzw",
               bg = "white", res = 300, type = "windows")

par(mfrow=c(1,2))
simplehist(kk[kk<4], xlab="Number of diseases affecting a farm per year", cex.lab=8/12,
           yaxt="n", xaxt="n")
axis(2, at=c(0, 2e4, 4e4, 6e4), labels=c(0, 2e4, 4e4, 6e4),
     cex.axis=8/12)
axis(1, at=c(0, 1, 2, 3), labels=c(0, 1, 2, 3),
     cex.axis=8/12)
simplehist(d$No_outbreaks, xlab="Number of diseases affecting a farm per year", 
           cex.lab=8/12, cex.axis=8/12)

dev.off()

length(kk[kk>3])/length(kk)
mean(kk)
var(as.vector(kk))
mean(d$No_outbreaks)
var(d$No_outbreaks)


# SImulate counts for poisson reg
# Samples per lambda:
k <- 1

pp <- matrix(rep(NA,8000*k*17), ncol=17)
for (i in 1:17){ 
  for (j in 1:8000){
  pp[((j-1)*k+1):(j*k),i] <- rpois(k,exp(lambda[j,i]))
  }
}

par(mfrow=c(1,2))
simplehist(pp[pp<4])
simplehist(d$No_outbreaks)
length(pp[pp>3])/length(pp)
mean(pp)
var(as.vector(pp))
length(pp[pp<1])/length(pp)
length(kk[kk<1])/length(kk)
length(d$No_outbreaks[d$No_outbreaks<1])/length(d$No_outbreaks)

# Further model checking
# Plotting
# No outbreaks

rate.med <- apply(kk, 2, FUN=median)
rate.90 <- apply(kk, 2,FUN=PI ,prob=0.90)

plot(rate.med ~ seq(1:nrow(d)), xlab="Site", ylab="Number of diseases",
     ylim=c(0, max(rate.90)), xaxt = "n")
axis(side = 1, at = seq(1:nrow(d)), lwd=0.5,
     labels = seq(1:nrow(d)))


for ( i in 1:nrow(d) ) {
  ci <- rate.90[,i]
  x <- seq(1:nrow(d))[i]
  lines( c(x,x) , ci)
  points(c(x,x) , ci, cex=0.7, pch=3)
}

points(seq(1:nrow(d)), d$No_outbreaks, cex=1,pch=20, col=rangi2)





# Inferences from top model: simulate to make inference and plot
## Comparing high and low contact for fixed high and low biosec
a<- as.vector(Intercept)
bs <- as.vector(Score)
bi <- as.vector(NoProviders)
bsi <- as.vector(Interaction)
lambda <- as.vector(lambda)

tiff(filename = "Fig3.tiff",
               width = 7.5, height = 4, units = "in", pointsize = 12,
               compression = "lzw",
               bg = "white", res = 325, type = "windows")

par(mfrow=c(2,2), mai=c(0.5, 0.3, 0.1, 0.3))
# Low biosec
dis_high_c <- exp(a + bi*max(d$windeg_c) + (bsi*max(d$windeg_c) + bs)*min(d$Score_surv100_c))
dis_low_c <- exp(a + bi*min(d$windeg_c) + (bsi*min(d$windeg_c) + bs)*min(d$Score_surv100_c))

# Difference in number of diseases
dis_diff <- dis_high_c - dis_low_c
sum(dis_diff > 0)/length(dis_diff)
plot(density(dis_diff), xlab="Difference high vs low No of suppliers | low biosec", main=""
     , xlim=c(-4, 4), cex.lab=8/12, cex.axis=8/12,  mgp = c(2, 1, 0))
shade(density(dis_diff), PI(dis_diff, prob=0.95))
abline(v=median(dis_diff), lty=4)
#a.1 <- locator()
text(a.1, "a)")
# High biosec
dis_high_c <- exp(a + bi*max(d$windeg_c) + (bsi*max(d$windeg_c) + bs)*max(d$Score_surv100_c))
dis_low_c <- exp(a + bi*min(d$windeg_c) + (bsi*min(d$windeg_c) + bs)*max(d$Score_surv100_c))

# Difference in number of diseases
dis_diff <- dis_high_c - dis_low_c
sum(dis_diff > 0)/length(dis_diff)
plot(density(dis_diff), xlab="Difference high vs low No of suppliers | high biosec", main=""
     ,xlim=c(-4, 4), cex.lab=8/12, cex.axis=8/12,  mgp = c(2, 1, 0))
shade(density(dis_diff), PI(dis_diff, prob=0.95))
abline(v=median(dis_diff), lty=4)
#b <- locator()
text(b, "b)")

## Comparing high and low biosec for fixed high and low contact rate

# High contact rate
dis_high_bs <- exp(a + bi*max(d$windeg_c) + (bsi*max(d$windeg_c) + bs)*max(d$Score_surv100_c))
dis_low_bs <- exp(a + bi*max(d$windeg_c) + (bsi*max(d$windeg_c) + bs)*min(d$Score_surv100_c))

# Difference in number of diseases
dis_diff <- dis_low_bs - dis_high_bs
sum(dis_diff > 0)/length(dis_diff)
plot(density(dis_diff), xlab="Difference low vs high biosecurity | high No of suppliers", main="",
     xlim=c(-4, 4), cex.lab=8/12, cex.axis=8/12,  mgp = c(2, 1, 0))
shade(density(dis_diff), PI(dis_diff, prob=0.95))
abline(v=median(dis_diff), lty=4)
#c <- locator()
text(c, "c)")

# Low contact rate
dis_high_bs <- exp(a + bi*min(d$windeg_c) + (bsi*min(d$windeg_c) + bs)*max(d$Score_surv100_c))
dis_low_bs <- exp(a + bi*min(d$windeg_c) + (bsi*min(d$windeg_c) + bs)*min(d$Score_surv100_c))

# Difference in number of diseases
dis_diff <- dis_low_bs - dis_high_bs
sum(dis_diff > 0)/length(dis_diff)
plot(density(dis_diff), xlab="Difference low vs high biosecurity | low No of suppliers", main="",
     xlim=c(-4, 4), cex.lab=8/12, cex.axis=8/12,  mgp = c(2, 1, 0))
shade(density(dis_diff), PI(dis_diff, prob=0.95))
abline(v=median(dis_diff), lty=4)
#d.1 <- locator()
text(d.1, "d)")

dev.off()

# Doing counterfactual predictions


# compute trend for high and low contact rate site
biosec.seq <- seq(from=min(d$Score_surv100_c)-.032, to=max(d$Score_surv100_c)+.032, length.out=80)


theta.pred.h <- matrix(rep(NA, 80*8000),ncol=80)


for (i in 1:length(biosec.seq)){
  theta.pred.h[,i] <- exp(a + bi*max(d$windeg_c) + (bsi*max(d$windeg_c) + bs)*biosec.seq[i])
}
  
theta.med.h <- apply( theta.pred.h, 2 , median )
theta.66.h <- apply( theta.pred.h, 2 , PI, prob=0.66 )
theta.90.h <- apply( theta.pred.h, 2 , PI, prob=0.95 )

theta.pred.l <- matrix(rep(NA, 80*8000),ncol=80)

for (i in 1:length(biosec.seq)){
  theta.pred.l[,i] <- exp(a + bi*min(d$windeg_c) + (bsi*min(d$windeg_c) + bs)*biosec.seq[i])
}

theta.med.l <- apply( theta.pred.l, 2 , median )
theta.66.l <- apply( theta.pred.l, 2 , PI, prob=0.66 )
theta.90.l <- apply( theta.pred.l, 2 , PI, prob=0.95 )

# Plotting
tiff(filename = "Fig2.tiff",
               width = 5.2, height = 5.2, units = "in", pointsize = 12,
               compression = "lzw",
               bg = "white", res = 350, type = "windows")

par(mfrow=c(1,1))
pch <- ifelse( d$windeg_c > 0 , 16 , 1 ) 
set.seed(0)
plot( jitter(d$Score_surv100_c, factor=15) , jitter(d$No_outbreaks,
      factor=0.0) , col=rangi2 , pch=pch , xlab="Biosecurity score",
      ylab="No of different diseases", ylim=c(0,3.0), xlim=c(-0.85,0.65),
      xaxt = "n",cex.lab=8/12, cex.axis=8/12,  mgp = c(2, 1, 0))
axis(side = 1, cex.axis=8/12, at = seq(from=min(d$Score_surv100_c), to=max(d$Score_surv100_c), length.out=11), lwd=0.5,
     labels =round(seq(from=79.5, to=88.6, length.out=11), digits=0))

lines( biosec.seq , theta.med.h , col=rangi2 )
shade( theta.90.h , biosec.seq , col=col.alpha("pink",0.2) )
shade( theta.66.h , biosec.seq , col=col.alpha(rangi2,0.2) )

lines( biosec.seq , theta.med.l , lty=2 )
shade( theta.90.l , biosec.seq , col=col.alpha("green",0.2) )
shade( theta.66.l , biosec.seq , col=col.alpha("black",0.2) )

dev.off()

# Does it add up to one? RMK calculates normalizing constant, so no need to normalize

pLGP(0:4,theta[1,1],lambda[1])
qLGP(c(0.01, 0.25, 0.5, 0.75, 1),theta[1,1],lambda[1])
dLGP(0:4,theta[1,1],lambda[1])
LGP.findmax(theta[1,1],lambda[1])
LGP.get.nc(theta[1,1],lambda[1])
LGPMVP(mu=mean(d$No_outbreaks),sigma2=var(d$No_outbreaks))


#########################################################################################
######################## Mortality analysis #############################################
#########################################################################################


# Doing mortality analysis
library(rethinking)
library(doBy)
mort <- read.csv("Mort_All_2012_2014.csv")

#offset and variables
mort$Closing.Count[is.na(mort$Closing.Count)] <- mort$Opening.Count[which.max(is.na(mort$Closing.Count))]
logpop <- log((mort$Opening.Count+ mort$Closing.Count)/2)
weight <- mort$Opening.Weight
mortality <- mort$Mortality.number
biosec <- scale(mort$Biosec, scale=sd(mort$Biosec, na.rm=T)*2)
windeg_c <- ifelse(mort$Suppliers>1, 1,0)
windeg_c <- windeg_c-mean(windeg_c)
inter <- biosec*windeg_c
site <- mort$Site
cage <- mort$Cage
yc <- mort$Yearclass
week <- ifelse(mort$Year==2012, mort$Week, ifelse(mort$Year==2013, 52 + mort$Week, 
                                                  52 + 52 + mort$Week))

mort2 <- data.frame(list(logpop=logpop, weight=weight, mort = mortality, biosec=biosec, windeg_c=windeg_c, 
                         inter=inter, site=site, cage=cage, yc=yc, week=week))

mort3 <- summaryBy(mort ~ site + cage + yc + week, FUN=sum, data=mort2)
logpop3 <- summaryBy(logpop + windeg_c + biosec + inter + weight ~ site + cage + yc + week, FUN=max, data=mort2)
mort3$logpop <- logpop3$logpop.max
mort3$weight <- logpop3$weight.max
mort3$windeg_c <- logpop3$windeg_c.max
mort3$biosec <- logpop3$biosec.max
mort3$inter <- logpop3$inter.max 

colnames(mort3) <- c("site", "cage", "yc", "week", "mort", "logpop",
                     "weight", "windeg_c", "biosec", "inter")



# Distance matrix
d.cent$sample.no = as.integer(as.factor(d.cent$Week))

# Creatig a pairwise difference in time matrix (time distance b/w obs)
dist.u <- dist(unique((mort3$week), method = "maximum", diag = T, upper = T, p = 2))
dist.m.u <- as.matrix(dist.u)
dist.m.u <- scale(dist.u, center=F, scale=rep(max(dist.u), ncol(dist.m.u)))




data.mort=list(logpop=mort3$logpop, mort = mort3$mort, site=mort3$site, cage=mort3$cage,
               yc=mort3$yc, week=mort3$week, N=nrow(mort3), windeg_c=mort3$windeg_c,
               biosec=mort3$biosec, inter=mort3$inter, weight=mort3$weight, Dmat=dist.m.u,
               N_site=length(unique(mort3$site)), N_cage=length(unique(mort3$cage)),
               N_yc=length(unique(mort3$yc)), N_week=length(unique(mort3$week)))



#No measurement error

m1.0 <- map2stan(
  alist(mort ~ dgampois( lambda , theta),
        log(lambda) <- logpop + a + a_site[site] + a_cage[cage],
        a ~ dnorm(0,10),
        theta ~ dcauchy(0,1),
        a_site[site] ~dnorm(0, sigma_site),
        a_cage[cage] ~dnorm(0, sigma_cage),
        c(sigma_site, sigma_cage) ~ dcauchy(0,1)
  ), data=data.mort,start=list(sigma_site=1.0, sigma_cage=1.0),
  iter=100 , warmup=50 , chains=1 ,cores=1,
  control = list(adapt_delta = 0.80, stepsize=1))

m1.0 <- map2stan(m1.0, data=data.mort,
                 iter=2000 , warmup=1000 , chains=2 ,cores=2,
                 control = list(adapt_delta = 0.80, stepsize=1))


precis(m1.0, depth=2)
plot(m1.0)
WAIC(m1.0)
postcheck(m1.0)



# Non-cenetered parameterization of varying effects
data.mort$site<- m1.0@data$site
data.mort$cage<- m1.0@data$cage
data.mort$yc<- m1.0@data$yc
data.mort$week<- m1.0@data$week

stanmort1 <- 
"data{
    int<lower=1> N;
int<lower=1> N_site;
int<lower=1> N_cage;
int<lower=1> N_yc;
int<lower=1> N_week;
int mort[N];
int site[N];
int cage[N];
int yc[N];
int week[N];
real logpop[N];
}
parameters{
real a;
real<lower=0> theta;
vector[N_site] a_site_raw;
vector[N_cage] a_cage_raw;
vector[N_yc] a_yc_raw;
vector[N_week] a_week_raw;
real<lower=0> sigma_site;
real<lower=0> sigma_cage;
real<lower=0> sigma_yc;
real<lower=0> sigma_week;
}
transformed parameters{
vector[N_site] a_site;
vector[N_cage] a_cage;
vector[N_yc] a_yc;
vector[N_week] a_week;

a_site = a_site_raw*sigma_site;
a_cage = a_cage_raw*sigma_cage;
a_yc = a_yc_raw*sigma_yc;
a_week = a_week_raw*sigma_week;
}

model{
vector[N] lambda;
sigma_week ~ cauchy( 0 , 2 );
sigma_yc ~ cauchy( 0 , 2 );
sigma_cage ~ cauchy( 0 , 2 );
sigma_site ~ cauchy( 0 , 2 );
a_week_raw ~ normal(0,1);
a_yc_raw ~ normal( 0 , 1 );
a_cage_raw ~ normal( 0 , 1 );
a_site_raw ~ normal( 0 , 1 );
theta ~ cauchy( 0 , 2 );
a ~ normal( 0 , 10 );
for ( i in 1:N ) {
lambda[i] = logpop[i] + a + a_site[site[i]] + a_cage[cage[i]] + a_yc[yc[i]] + a_week[week[i]];
lambda[i] = exp(lambda[i]);
}
mort ~ neg_binomial_2( lambda , theta );
}
generated quantities{
vector[N] lambda;
vector[N] log_lik;

for ( i in 1:N ) {
lambda[i] = logpop[i] + a + a_site[site[i]] + a_cage[cage[i]] + a_yc[yc[i]] + a_week[week[i]];
lambda[i] = exp(lambda[i]);
log_lik[i] = neg_binomial_2_lpmf( mort[i] | lambda[i] , theta );

}
}
"


m1.1 <- stan(model_code = stanmort1, iter=2000, chains=2, cores=2,
            warmup=1000, control = list(adapt_delta = 0.95),
            data=data.mort)

print(m1.1, probs=c(0.025, 0.975), pars=c("a_site_raw", "a_cage_raw", "a_yc_raw", "a_week_raw", 
                                          "log_lik", "lambda"), include=F)


traceplot(m1.1, pars=c("a", "a_site"))
traceplot(m1.1, pars= paste(rep("a_cage", 25), rep("[", 25),seq(from=1, to=25), rep("]", 25), sep = ""))
traceplot(m1.1, pars= paste(rep("a_cage", 25), rep("[", 25),seq(from=26, to=50), rep("]", 25), sep = ""))
traceplot(m1.1, pars= paste(rep("a_cage", 25), rep("[", 25),seq(from=51, to=75), rep("]", 25), sep = ""))
traceplot(m1.1, pars= paste(rep("a_cage", 25), rep("[", 25),seq(from=76, to=100), rep("]", 25), sep = ""))
traceplot(m1.1, pars= paste(rep("a_cage", 25), rep("[", 25),seq(from=101, to=125), rep("]", 25), sep = ""))
traceplot(m1.1, pars= paste(rep("a_cage", 25), rep("[", 25),seq(from=126, to=150), rep("]", 25), sep = ""))
traceplot(m1.1, pars= paste(rep("a_cage", 25), rep("[", 25),seq(from=151, to=171), rep("]", 25), sep = ""))
traceplot(m1.1, pars=c("a_yc"))
traceplot(m1.1, pars= paste(rep("a_week", 25), rep("[", 25),seq(from=1, to=25), rep("]", 25), sep = ""))
traceplot(m1.1, pars= paste(rep("a_week", 25), rep("[", 25),seq(from=26, to=50), rep("]", 25), sep = ""))
traceplot(m1.1, pars= paste(rep("a_week", 25), rep("[", 25),seq(from=51, to=75), rep("]", 25), sep = ""))
traceplot(m1.1, pars= paste(rep("a_week", 25), rep("[", 25),seq(from=76, to=100), rep("]", 25), sep = ""))
traceplot(m1.1, pars= paste(rep("a_week", 25), rep("[", 25),seq(from=101, to=125), rep("]", 25), sep = ""))
traceplot(m1.1, pars= paste(rep("a_week", 25), rep("[", 25),seq(from=126, to=150), rep("]", 25), sep = ""))
traceplot(m1.1, pars= paste(rep("a_week", 25), rep("[", 5),seq(from=151, to=155), rep("]", 5), sep = ""))
traceplot(m1.1, pars=c("sigma_site","sigma_cage", "sigma_yc", "sigma_week"))


library(loo)
loglikmort <- extract_log_lik(m1.1)
(mort_loo <- loo(loglikmort))


pairs(mkk, pars=c("bagd", "bdes", "bcyst", "bprv", "bpox", "bten", "bagd_des", "bagd_cyst", "bagd_prv", "bagd_pox", 
                  "bagd_ten","btemp", "btime_w", "btemp_time"))




stan_plot(mkk, pars=c("bagd", "bdes", "bcyst", "bprv", "bpox", "bten", "bagd_des", "bagd_cyst", "bagd_prv", "bagd_pox", 
                      "bagd_ten","btemp", "btime_w", "btemp_time"))
