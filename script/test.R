library(tidyverse)
loadNamespace("cowplot")
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
longer.selected.vars <- selected.vars %>%  arrange(X1) %>% 
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


dbb.25.data <- read_tsv("dbb.tsv")


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



    


before.data <- read_tsv("before.tsv")
after.data <- read_tsv("after.tsv")
em.data <- read_tsv("em.tsv")


bind_rows(after.data %>% select(purity) %>% mutate(type="after"),
          before.data %>% select(purity) %>% mutate(type="before")) %>%
    ggplot() + geom_histogram(aes(x=purity),bins=100) + facet_grid(type~.)

pvalues <- read_tsv("pvalues.tsv")

pval.unit <- full_join(pvalues, before.data, by="unit")


tig.data <- read_tsv("dbb.contig")
unit.data <- read_tsv("dbb.unit")

tig.data <- read_tsv("test.contig")

aln.tig006 <- read_csv("tig_006.tsv")
dbb25.contact <- read_tsv("dbb.25.contact.tsv")

lc.data <- read_tsv("./check.tsv", col_names = c("RESULT","seed","TYPE","rand","time","acc"))
lc.data <- read_tsv("./check_sr.tsv", col_names = c("RESULT","seed","TYPE","rand","time","acc"))

homopolymer_count <- read_tsv("dbb.tsv")
homopolymer_count %>% filter(Len > 10) %>%
    ggplot() + geom_point(aes(x=Len,y=Count,color=Base))

homopolymer_count <- bind_rows(read_tsv("dbb.tsv") %>% mutate(ds = "dbb"),
                               read_tsv("hg002.tsv") %>% mutate(ds = "hg002"),
                               read_tsv("lerkyo.tsv") %>% mutate(ds="lerkyo"))

homopolymer_count %>% filter(Len > 50) %>%
    ggplot() + geom_line(aes(x=Len,y=Count,color=Base)) + facet_grid(ds~.)


graph.features <- read_tsv("dip-A-CS.dump.tsv")
graph.features <- graph.features %>% mutate(skewness = hap1/(hap1+hap2)-1/2)

cov.unit.data <- read_tsv("units.tsv")
cov.flip.data <- read_tsv("flip.tsv")
cov.flip.unit <- full_join(cov.flip.data,
                           cov.unit.data, 
                           by = c("unit","cluster"))
consis.data <- read_tsv("consis.tsv") %>% mutate(idx = iter * 2 + loop)

consis.data  %>% filter(iter > 5) %>% ggplot() + geom_line(aes(x=iter,y=consis, color=factor(loop)))

cov.flip.unit %>% ggplot(aes(x=hap1, y=hap2)) + geom_point() + facet_wrap(vars(phaseblock))

cov.flip.unit %>% ggplot(aes(x=hap1, y=hap2, color = factor(phaseblock))) + geom_point()


cov.flip.unit %>% group_by(phaseblock) %>% summarize(n=n(), hap1=sum(hap1), hap2=sum(hap2))

cov.flip.unit %>% filter(copynum >1) %>% group_by(unit) %>%
    summarize(flip=sum(flip), copy_num=min(copynum), purity = sum(purity)/n()) %>%
    filter(flip < 20 & purity > 0.8 & purity < 1)


cov.flip.unit %>% filter(copynum >1) %>% group_by(unit) %>% filter(n() > 1) %>% 
    summarize(flip=sum(flip), copy_num=min(copynum), purity = sum(purity)/n()) %>%
    ggplot() + geom_point(aes(x=flip, y = purity))


cov.flip.unit %>% filter(copynum >1) %>% group_by(unit) %>% filter(n() > 1) %>% 
    summarize(flip=sum(flip), copy_num=min(copynum), purity = sum(purity)/n()) %>%
    ggplot() + geom_point(aes(x=flip, y = purity))


cov.flip.unit.summary <- cov.flip.unit %>% group_by(unit) %>% 
    summarize(flip=sum(flip, na.rm = TRUE),
              copy_num=min(copynum, na.rm = TRUE),
              purity = mean(purity, na.rm = TRUE),
              score = mean(score, na.rm = TRUE))
cov.flip.unit.summary %>%
    filter(score < 500 & flip < 75) %>% 
    ggplot() + geom_point(aes(x=flip, y = score, color = purity)) + scale_color_gradient2(high="blue", low="orange", midpoint = 0.8) 


cov.flip.unit.summary %>%  ggplot() + geom_histogram(aes(x=flip))
cov.flip.unit.summary %>% filter(score < 500) %>%  ggplot() + geom_histogram(aes(x=score), bins=100)



error.rate.data <- read_tsv("./result/sim.error.tsv")%>% mutate(total = mism + ins + del) 

error.rate.data %>% ggplot() + geom_histogram(aes(x=mism+ins+del))

error.rate.data %>% summarize(across(everything(), sd))
puerror.rate.data %>% summarize(across(everything(), mean))

error.rate.data %>% mutate(error = mism + ins + del) %>% select(-mism,-ins,-del) %>%
    group_by(unitid) %>% summarize(mean = mean(error), sd =sd(error)) %>%
    ggplot() + geom_point(aes(x=mean, y = sd/mean))

error.rate.data %>% mutate(error = mism + ins + del) %>% select(-mism,-ins,-del) %>%
    group_by(readid) %>% summarize(mean = mean(error), sd =sd(error)) %>%
    ggplot() + geom_point(aes(x=mean, y = sd/mean))


error.rate.data <- read_tsv("./result/real.error.tsv")

error.rate.data %>% ggplot() + geom_histogram(aes(x=mism+ins+del))


error.rate.data %>% group_by(unit,cluster) %>% summarize(total=mean(mism+ins+del))

error.rate.data %>% mutate(total = mism + ins + del) %>% summarize(across(everything(), sd))
error.rate.data %>% mutate(total = mism + ins + del) %>% summarize(across(everything(), mean))


error.rate.data %>% mutate(error = mism + ins + del) %>% select(-mism,-ins,-del) %>%
    group_by(unitid) %>% summarize(mean = mean(error), sd =sd(error)) %>%
    ggplot() + geom_point(aes(x=mean, y = sd/mean))

error.rate.data %>% mutate(error = mism + ins + del) %>% select(-mism,-ins,-del) %>%
    group_by(readid) %>% summarize(mean = mean(error), sd =sd(error)) %>%
    ggplot() + geom_point(aes(x=mean, y = sd/mean))



datasets<- map(list.files("./temp/", pattern = ".tsv"), ~ list(name=., data=read_tsv(paste0("./temp/",.))))


summarize_dataset <- function(df){
    df %>% group_by(unit) %>%
        filter(1 < n()) %>% ungroup() %>% summarize(n=n(), purity=mean(purity))
}

summarize_to_df <- function(ls){
    summarize_dataset(ls$data) %>% mutate(id=gsub("COX_PGF_ONT_30x_aligned_reads.","",ls$name))
}

datasets %>% map_dfr(summarize_to_df)

correction.comparison <- full_join(
    datasets[[6]]$data %>% group_by(unit) %>% summarize(purity=mean(purity)) %>% rename(after=purity),
    datasets[[2]]$data %>% group_by(unit) %>% summarize(purity=mean(purity)) %>% rename(before=purity))


error.pre <- read_tsv("./result/error.pre.tsv")
error.post <- read_tsv("./result/error.post.tsv")
