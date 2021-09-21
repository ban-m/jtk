library(tidyverse)
xlogx <- function(x) ifelse(abs(x)<0.0001, 0, x*log(x))

corel_model_aic <- function(xs,ys){
    x1_true_prob <- sum(xs) / length(xs)
    error_prob <- sum(xs == ys) / length(xs)
    error_prob <- max(error_prob, 1-error_prob)
    partlk <- function(p) xlogx(p)+xlogx(1-p)
    -2 * (partlk(x1_true_prob) + partlk(error_prob)) * length(xs) - 2 * 2
}


nest_model_aic <- function(xs,ys){
    x1_true_prob <- sum(xs) / length(xs)
    x1_t_error_prob <- sum((xs == ys) & xs) / sum(xs)
    x1_t_error_prob <- max(x1_t_error_prob, 1-x1_t_error_prob)
    x1_f_error_prob <- sum((xs == ys) & !xs) / sum(!xs)
    x1_f_error_prob <- max(x1_f_error_prob, 1- x1_f_error_prob)
    partlk <- function(p) xlogx(p)+xlogx(1-p)
    lk <- partlk(x1_true_prob) * length(xs) + partlk(x1_t_error_prob) * sum(xs) + partlk(x1_f_error_prob) * sum(!xs)
    -2 * lk - 2 * 3
}

full_model_aic <- function(xs,ys){
    tt <- sum((xs == ys) & xs)
    ff <- sum(xs==ys) - tt
    tf <- sum((xs != ys) & xs)
    ft <- sum(xs != ys) - tf
    totlen <- length(xs)
    lk <- (xlogx(tt/totlen) + xlogx(tf/totlen) + xlogx(ft/totlen) + xlogx(ff/totlen)) * totlen
    -2 * lk - 2 * 3
}

mat <- read_tsv("pmat.tsv")
selected.vars <- read_tsv("selected.tsv", col_names = FALSE)
selected.vars <- read_tsv("selected.2.tsv", col_names = FALSE)


len <- length(selected.vars$X1)
longer.selected.vars <- selected.vars %>% # arrange(X1) %>% 
    mutate(ReadID=1:len) %>% pivot_longer(cols = !ReadID & !X1, names_to = "Position", values_to = "LK") %>% 
    mutate(Position = as.integer(str_sub(Position,start=2)))


all.read.profiles <- read_tsv("pmat.tsv")
selected.profiles <- read_tsv("pos.tsv",col_names=c("Pos","Type","LK"))

into.base.type <- function(x){(c("A","C","G","T"))[x %% 4 + 1 ]}
into.edit.type <- function(x){(c("Sub","Ins","Del"))[x %/% 4 + 1]}


top.2.variants <- selected.profiles %>%
    mutate(EditOp = into.edit.type(Type), BaseType = into.base.type(Type)) %>%
    tail(n=2) %>% transpose()


temp <- all.read.profiles %>% filter(
(Position == 1418 & EditType == "Sub" & Base == "A") |
(Position == 320 & EditType == "Sub" & Base == "G")
)  %>% select(ReadID, Position, LK) %>%  pivot_wider(names_from = Position, values_from = LK) %>%
    rename("Var1"= `320`, "Var2" = `1418`)


temp <- all.read.profiles %>% filter(
(Position == 1418 & EditType == "Sub" & Base == "A") |
(Position == 1555 & EditType == "Ins" & Base == "G")
)  %>% select(ReadID, Position, LK) %>%  pivot_wider(names_from = Position, values_from = LK) %>%
    rename("Var1"= `1555`, "Var2" = `1418`)


temp <- all.read.profiles %>% filter(
(Position == 1418 & EditType == "Sub" & Base == "A") |
(Position == 55 & EditType == "Sub" & Base == "G")
)  %>% select(ReadID, Position, LK) %>%  pivot_wider(names_from = Position, values_from = LK) %>%
    rename("Var1"= `55`, "Var2" = `1418`)


error_rate <- 0.05
size = 80
coverage = 27
poismax <- function(n){
    max(1:10 %>% (function(c){dpois(x = n, lambda = coverage * c)}))
}

error_thresholds <- tibble(n=0:40) %>% mutate(binomprob = dbinom(x=n,size = size, prob=2*error_rate*(1-error_rate)),
                          poisprob = sapply(X=n, FUN=poismax))
    
    
IsTheSameBipart <- function(error.num,length, coverage,cl.error){
    pois.prob <- max(1:10 %>% (function(c){dpois(x = error.num, lambda = coverage * c)}))
    binom.prob <- dbinom(x=error.num, prob = 2 * cl.error * (1-cl.error), size = length)
    print(pois.prob)
    print(binom.prob)
    pois.prob < binom.prob
}


selected.vars <- read_tsv("selected.tsv", col_names = FALSE)
len <- length(selected.vars$X1)
longer.selected.vars <- selected.vars %>%  arrange(X1) %>% select(-X1) %>%  
    mutate(ReadID=1:len) %>% pivot_longer(cols = !ReadID ,names_to = "Position", values_to = "LK") %>% 
    mutate(Position = as.integer(str_sub(Position,start=2)))
g<- longer.selected.vars %>% ggplot() + geom_raster(aes(x=Position, y = ReadID, fill=LK)) +
    scale_fill_gradient2(low="blue", high="orange")




longer.selected.vars <- selected.vars %>%
    arrange(X2) %>% 
    mutate(ReadID=1:len) %>% pivot_longer(cols = !ReadID & !X1, names_to = "Position", values_to = "LK") %>% 
    mutate(Position = as.integer(str_sub(Position,start=2)))
g<- longer.selected.vars %>% ggplot() + geom_raster(aes(x=Position, y = ReadID, fill=LK)) + scale_fill_gradient2(low="blue", high="orange")


selected.vars %>% group_by(X1) %>%
    summarise(across(everything(), sum))

total.num <- 80
true.prob <- 0.25
bind_rows(
    tibble(x=1:40) %>% mutate(y = dpois(x, lambda = total.num * true.prob, log = TRUE), type = "Poisson"),
    tibble(x=1:40) %>% mutate(y = dbinom(x, size = total.num, p=true.prob, log=TRUE), type = "Binom")
)%>% ggplot() + geom_line(aes(x=x,y=y,color=type))


cosine <- function(xs,ys){
    sum(xs * ys) / sqrt(sum(xs^2)) / sqrt(sum(ys^2))
}



clustering.num <- read_tsv("clustering.tsv", col_names = c("ID","Time","Before","After"))

assignments <- read_tsv("assignments.tsv", col_names = FALSE)

variants <- read_tsv("vars.tsv", col_names=c("read", "position", "lk"))
variants <- read_tsv("vars.normal.tsv", col_names=c("read", "position", "lk"))

variants <- variants %>%
    mutate(location = position %/% 9, edit_type = (position %% 9) %/% 4, base = (position %%9 ) %% 4)


variants %>% ggplot() + geom_raster(aes(y=read, x=location, fill=lk)) + scale_fill_gradient2(low="blue", high="orange") + facet_grid(edit_type~base)


variants %>%
    filter(900 < location & location < 1100) %>%
    ## filter(base == 0 | base == 3) %>%
    mutate(lk = ifelse(lk>0,lk,rep(0,length(lk)))) %>% 
    ggplot() + geom_raster(aes(y=read, x=location, fill=lk)) + scale_fill_gradient2(low="blue", high="orange") + facet_grid(edit_type~base)

normalize <- function(lk){
    mean <- mean(lk)
    sd <- sd(lk)
    ifelse(rep(is.na(sd),length(lk)), rep(0,length(lk)), (lk-mean)/sd)
}

normalized.variants <- variants %>% group_by(position) %>%
    nest() %>% 
    mutate(data=map(data,~ mutate(.x, lk = (lk - mean(lk))/ ))) %>%
    unnest(cols = c(data))

normalized.variants %>%
    filter(990 < location & location < 1050) %>% 
    ggplot() + geom_raster(aes(y=read, x=location, fill=lk)) + scale_fill_gradient2(low="blue", high="orange") + facet_grid(edit_type~base)

normalized.variants.sd <- variants %>% group_by(position) %>%
    nest() %>% 
    mutate(data=map(data,~ mutate(.x, lk = normalize(lk)))) %>%
    unnest(cols = c(data))

normalized.variants.sd %>%
    ggplot() + geom_raster(aes(y=read, x=location, fill=lk)) + scale_fill_gradient2(low="blue", high="orange") + facet_grid(edit_type~base)

features.vector <- normalized.variants.sd %>% group_by(position) %>% nest() %>%
    mutate(data=map(data,~summarize(.x, location=location[1], edit_type=edit_type[1], base = base[1], vars = sum(lk[lk>0]), num = sum(lk>0)))) %>%
    unnest(cols=c(data))



combined.data <- corel.data %>% group_by(ID1) %>% nest() %>%
    mutate(data=map(data,~summarize(.x,PValue=min(PValue)))) %>%
    unnest(cols=c(data)) %>%
    rename("unit"=ID1) %>% 
    full_join(unit.pre, by="unit")

features %>% filter(type=="Cp" & pos < 500)  %>% ggplot() + geom_raster(aes(x=pos,y=readid,fill=lk)) + scale_fill_gradient2(low="blue",high = "orange") + facet_wrap(.~len) + labs(x="Position(bp)", y="ReadID", fill="Fold change(Score)")

edge.data <- read_tsv("dip-A-CS.edge")
edge.data <- read_tsv("dbb.edge")
edge.data <- read_tsv("hg002.edge")
g <- edge.data %>% filter(len < 6000) %>%  ggplot() + geom_histogram(aes(x=len))


tig.data <- read_tsv("sq.contig", col_names = c("id","copy","len","hap1","hap2")) %>%
    mutate(total=hap1+hap2, acc=pmax(hap1/total,hap2/total))


tig.data <- read_tsv("pre_em.tsv", col_names = c("id","copy","len","hap1","hap2")) %>%
    mutate(total=hap1+hap2, acc=pmax(hap1/total,hap2/total))


em.lk.data <- read_tsv("em_data.tsv")

em.lk.summary <- em.lk.data %>% group_by(num,span) %>% nest() %>%
    mutate(data=map(data,~summarize(.x,lk=mean(lk),null=mean(null)))) %>%
    unnest(cols=c(data))

em.lk.summary <- em.lk.data %>% mutate(diff=lk-null) %>% group_by(span,num) %>%
    nest() %>%
    mutate(data=map(data,~summarize(.x,  counts=sum(diff>0), diff=mean(diff), rand_index=mean(rand_index)))) %>%
    unnest(cols=c(data))



em.lk.summary.onlysuc <- em.lk.data %>% filter(clusternum==2) %>% group_by(span,num) %>%
    nest() %>%
    mutate(data=map(data,~summarize(.x,  rand_index=mean(rand_index)))) %>%
    unnest(cols=c(data))

em.lk.summary %>% ggplot() + geom_line(aes(x=num,y=diff)) + facet_wrap(.~span)

span.probs <- read_tsv("span.tsv")


var_pos <- read_csv("diffs.csv",col_names=FALSE)
var_pos %>% ggplot() + geom_point(aes(x=X1,y=X2,size=X3,color=X4))


selected.answers <- read_tsv("75.out", col_names = FALSE) %>% select(X1,X2) %>% rename(answer=X2)
selected.vars <- read_tsv("75.out", col_names = FALSE) 
len <- length(selected.vars$X1)

longer.selected.vars <- selected.vars %>%
    rename(ReadID=X1)%>% 
    arrange(X2) %>%
    pivot_longer(cols= !ReadID, names_to = "Position", values_to = "LK") 
g<- longer.selected.vars %>% ggplot() + geom_raster(aes(x=Position, y = ReadID, fill=LK)) +
    scale_fill_gradient2(low="blue", high="orange")


longer.selected.vars <- selected.vars %>%
    arrange(X2) %>%
    rename(ReadID=X1)%>%
    mutate(ReadID=1:len) %>% 
    pivot_longer(cols= !ReadID, names_to = "Position", values_to = "LK") 
g<- longer.selected.vars %>% ggplot() + geom_raster(aes(x=Position, y = ReadID, fill=LK)) +
    scale_fill_gradient2(low="blue", high="orange")

scores <- read_tsv("./score.tsv")


pre.em.data <- read_tsv("pre_em.unit", col_names = c("id","cluster","hap1","hap2")) %>%
    mutate(total=hap1+hap2, acc=pmax(hap1/total,hap2/total))

em.data <- read_tsv("span.unit", col_names = c("id","cluster","hap1","hap2")) %>%
    mutate(total=hap1+hap2, acc=pmax(hap1/total,hap2/total))

rand.data <- read_tsv("sq.unit", col_names = c("id","cluster","hap1","hap2")) %>%
    mutate(total=hap1+hap2, acc=pmax(hap1/total,hap2/total))

joind.data <- full_join(em.data %>% rename(old=acc) %>% select(id,cluster,old) ,
                        rand.data %>% rename(new=acc) %>% select(id,cluster,new), by =c("id","cluster"))

full.data <- joind.data %>% left_join(y=pre.em.data,by=c("id","cluster")) 

merged.data <- full_join(rand.data %>% select(id,cluster,acc) %>% rename(em=acc),
                         pre.em.data %>% select(id,cluster,acc),
                         by = c("id", "cluster"))


diplotig.data <- read_tsv("diplotig.unit", col_names = c("id","cluster","hap1","hap2")) %>%
    mutate(total=hap1+hap2, acc=pmax(hap1/total,hap2/total))


dbb.25.data <- read_tsv("dbb.units")


sampling <- function(coverage,n){
    coverage <- rpois(lambda = coverage, n = n)
    sample_one_unit <- function(coverage){
        hap1 <- rbinom(size=coverage,n=1,prob=0.5)
        hap2 <- coverage - hap1
        list(hap1=hap1,hap2=hap2)
    }
    coverage %>% lapply(FUN=sample_one_unit)
}

test.data <- sampling(30,4000) %>% map_dfr(~tibble(hap1=.x$hap1, hap2=.x$hap2))

### Create (n,m)-[0,1]matrix with Pr{x_ij = 1} = p
sampling <- function(column, row, prob) {
    matrix(as.integer(rbernoulli(n=column * row, p=prob)),nrow=row)
}

max_gain <- function(column, row, prob){
    max(colSums(sampling(column,row,prob)))
}

parameters <- expand_grid(row=seq(from=20,to=90,by=10), column = seq(from = 1500,to=2000,by=100)) %>%
    mutate(count=300) %>% uncount(count)
probability <- 0.05
sim.result <- parameters %>% mutate(max = mapply(FUN = max_gain, column, row, MoreArgs = list(prob=probability)))


parameters <- expand_grid(row=seq(from=20,to=90,by=10), prob = seq(from=0.01,to=0.07,by=0.01), column=2000) %>% 
    mutate(count=300) %>% uncount(count)
sim.result <- parameters %>% mutate(max = mapply(FUN = max_gain, column, row, prob))


parameters <- expand_grid(row=seq(from=20,to=90,by=10), column = c(2000)) %>%
    mutate(count=300) %>% uncount(count)
probability <- 0.025
sim.result <- parameters %>% mutate(max = mapply(FUN = max_gain, column, row, MoreArgs = list(prob=probability)))



    

