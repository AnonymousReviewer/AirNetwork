
library(dplyr)
library(sp)
library(reshape2)
library(igraph)
library(Matrix)
library(ggplot2)

# data from https://www.bts.gov/browse-statistical-products-and-data/bts-publications/data-bank-28dm-t-100-domestic-market-data-us
raw <- as.tbl(read.csv('data/dd.db28dm.201701.201712.asc',sep='|',header=F,stringsAsFactors = F))
names(raw)<- c('year','month','oid','odigits','okey','oname','did','ddigits','dkey','dname','company','company_code',
               'type','odid','type2','passengers','freight','mail')

# group by origin/destination airports - checkodid -> statistically ok
summary(raw %>% group_by(oid,did) %>% summarize(ods = length(unique(odid))))

# sum passenger flows
flows <- raw %>% group_by(oid,did) %>% summarize(flow = sum(passengers))

# largest flow is LAX / SFO -> 1.7Mio passengers - OK see https://en.wikipedia.org/wiki/Los_Angeles_International_Airport#Passenger
flows[flows$flow == max(flows$flow),]
flows[flows$oid == "SFO"&flows$did=="LAX",]

# filter on flows with more than 10000 passengers
# -> done after aggregating both directions
#flows <- flows[flows$flow>10000,]
#length(unique(c(flows$oid,flows$did)))

# airport data from FAA
# https://ais-faa.opendata.arcgis.com/datasets/e747ab91a11045e8b3f8a3efd093d3b5_0/data
airports <- as.tbl(read.csv('data/Airports.csv',stringsAsFactors = F))
airports <- airports[airports$IDENT%in%flows$oid|airports$IDENT%in%flows$did,]

# construct the graph
nodes = airports[,c("IDENT","X","Y")]
nodenames = unique(c(flows$oid,flows$did))
missing = nodenames[!nodenames%in%nodes$IDENT]

g <- graph_from_data_frame(flows[!flows$oid%in%missing&!flows$did%in%missing,],directed = T,vertices = nodes)
E(g)$weight <- E(g)$flow
A = as_adjacency_matrix(g,attr = 'weight')
g<- graph_from_adjacency_matrix(A+t(A),weighted = T,mode = 'undirected')

g <- subgraph.edges(g,which(E(g)$weight>10000))

# great circle distances between nodes
#dmat = spDists(as.matrix(airports[,c("X","Y")]),longlat = T)
#rownames(dmat)<-airports$IDENT;colnames(dmat)<-airports$IDENT
#dists <- melt(dmat);names(dists)<-c("o","d","distance")
#flows <- left_join(flows,dists,by=c("oid"="o","did"="d"))

# directly distance between nodes
rownames(airports)<-airports$IDENT
V(g)$x <- unlist(airports[V(g)$name,"X"])
V(g)$y <- unlist(airports[V(g)$name,"Y"])

dmat = spDists(as.matrix(airports[V(g)$name,c("X","Y")]),longlat = T)
rownames(dmat)<-V(g)$name;colnames(dmat)<-V(g)$name
dists <- melt(dmat);names(dists)<-c("o","d","distance")
rownames(dists)<-paste0(dists$o,dists$d)

E(g)$distance <- dists[paste0(head_of(g,E(g))$name,tail_of(g,E(g))$name),"distance"]

# keep largest component only
g = induced_subgraph(g,which(components(g)$membership==1))

# average shortest path without weighting
mean(distances(g,weights = NA))

plot(g,vertex.label=NA,vertex.size=0.1)


###
# correlation between bw centrality and degree

# unweighted
cor.test(betweenness(g),degree(g))
# 0.79

# weighted
cor.test(betweenness(g,weights = NA),degree(g))
# 0.81
# -> very high, confirms intuition that in this kind of network, betweeness is not that relevant

###
# Same plot as Fig 4

E(g)$flow = E(g)$weight
E(g)$weight = E(g)$distance


# fixed variables for gbc
dists = distances(g)
nodepaths = lapply(V(g),function(v){shortest_paths(g,from=v,to=V(g))$vpath})

#'
#' nodesize should be named
#' rough implementation
computegbc <- function(beta,nodesize,sample=T){
  gbc = rep(0,length(nodesize));names(gbc)<-names(nodesize)
  egbc = rep(0,length(E(g)));names(egbc)<-paste0(head_of(g,E(g))$name,tail_of(g,E(g))$name)
  flows = matrix(rep(nodesize,length(nodesize)),nrow=length(nodesize),byrow = T) * matrix(rep(nodesize,length(nodesize)),nrow=length(nodesize),byrow = F) / dists^beta
  flows[is.infinite(flows)]=0
  # sum(flows[flows>mean(flows)])/sum(flows) = 0.94 : taking larger than the average makes 94% of total flow -> add filter
  meanflow = mean(c(flows))
  
  for(o in names(nodepaths)){
    show(o)
    for(d in 1:length(nodepaths[[o]])){
      path = nodepaths[[o]][[d]]
      oid = path[1]$name;did=path[length(path)]$name
      if(sample&flows[oid,did]>meanflow){
        #show(path)
        if(dists[oid,did]>0){
          flow = nodesize[oid]*nodesize[did]/(dists[oid,did]^beta)
          #for(n in 1:length(path)){gbc[path[n]$name]=gbc[path[n]$name]+flow}
          gbc[path$name]=gbc[path$name]+flow
          edgenames = paste0(path[1:(length(path)-1)]$name,path[2:length(path)]$name)
          revedgenames = paste0(path[2:length(path)]$name,path[1:(length(path)-1)]$name)
          #egbc[edgenames[edgenames%in%names(egbc)]] = egbc[edgenames[edgenames%in%names(egbc)]]+flow # too slow
          #egbc[revedgenames[revedgenames%in%names(egbc)]] = egbc[revedgenames[revedgenames%in%names(egbc)]]+flow
          newflow = egbc[edgenames]+flow;egbc[names(newflow)[!is.na(names(newflow))]] = newflow[!is.na(names(newflow))] # terrible
          newflow = egbc[revedgenames]+flow;egbc[names(newflow)[!is.na(names(newflow))]] = newflow[!is.na(names(newflow))]
        }
      }
    }
  }
  return(list(gbc=gbc,egbc=egbc))
}

weighteddegs = strength(g,weights = E(g)$flow)
beta=1

weightedgbc <- computegbc(beta,weighteddegs)

# correlation with degree for nodes
cor.test(weightedgbc$gbc,weighteddegs)

# correlation with passenger flow for edges
cor.test(weightedgbc$egbc,E(g)$flow)

dummypops = rep(1,length(weighteddegs));names(dummypops)<-names(weighteddegs)

bccor = cor.test(edge_betweenness(g),E(g)$flow)

betas = c();types = c();corrs = c();corrmins=c();corrmaxs=c()
for(beta in seq(0,2,0.2)){
  show(beta)
  weightedgbc <- computegbc(beta,weighteddegs)
  dummygbc <- computegbc(beta,dummypops)
  
  weightedcor = cor.test(weightedgbc$egbc,E(g)$flow)
  corrs=append(corrs,weightedcor$estimate);corrmins=append(corrmins,weightedcor$conf.int[1]);corrmaxs=append(corrmaxs,weightedcor$conf.int[2]);betas=append(betas,beta);types=append(types,"weighted")
  
  dummycor = cor.test(dummygbc$egbc,E(g)$flow)
  corrs=append(corrs,dummycor$estimate);corrmins=append(corrmins,dummycor$conf.int[1]);corrmaxs=append(corrmaxs,dummycor$conf.int[2]);betas=append(betas,beta);types=append(types,"dummy")
  
  corrs=append(corrs,bccor$estimate);corrmins=append(corrmins,bccor$conf.int[1]);corrmaxs=append(corrmaxs,bccor$conf.int[2]);betas=append(betas,beta);types=append(types,"bc")
  
}

res = data.frame(betas,types,corrs,corrmins,corrmaxs)

g=ggplot(res,aes(x=betas,y=corrs,color=types))
g+geom_point()+geom_line()+geom_errorbar(aes(ymin=corrmins,ymax=corrmaxs))+xlab(expression(beta))+
  ylab(expression(rho))
ggsave('correlations.png',width=20,height=15,units='cm')








