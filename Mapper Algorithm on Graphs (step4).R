Mapper=function(){
  library(foreach)
  library(doParallel)
  library(ggfortify)
  library(TDAmapper)
  library(igraph)
  library(umap)
  library(proxyC)
  library(cluster)
  library(bioDist)
  library(dbscan)
  num_int_seq = seq(10,80,5)
  percent_overlap_seq = seq(20,90,10)
  #num_bins_when_clustering_seq = seq(5,15,5)
  parameters=vector(mode="list",length=360)
  i=1
  for (interv in num_int_seq){
    for (per in percent_overlap_seq){
      #for (bins in num_bins_when_clustering_seq){
      vec=vector(length=2)
      names(vec)=c("interv","per")
      vec["interv"]=interv
      vec["per"]=per
      #vec["bins"]=bins
      parameters[[i]]=vec
      i=i+1
      #}
    }
  }
  
  prep=function(){
    dir.create("1")
    dir.create("2")
    dir.create("3")
    dir.create("4")
    dir.create("5")
    dir.create("6")
  }
  
  filter_comp <-function(d){
    return(umap(d)$layout)
  }
  
  prep()
  symbol=readline("Enter dataset's symbol:")
  for(l in 1:6){
    if(length(list.files(paste("./",l,sep="")))==120){
      next
    }
    print(paste("Current Timestep:",l))
    load(paste("./g",l,".RData",sep=""))
    dat=g
    rem(g)
    ppi_dist_matrix=dat
    emb_ppi=read.table(file=paste("./No_Feature/",symbol,"E",l,".csv",sep=""),header = TRUE,sep = ",",quote = "")
    rownames(emb_ppi)=emb_ppi[,1]
    emb_ppi=emb_ppi[,-1]
    emb_ppi=filter_comp(emb_ppi)
    All_Mapper_umap_p=vector(length=120,mode='list')
    totalCores=detectCores()
    clustery <- makeCluster(totalCores[1]-6)
    registerDoParallel(clustery)
    foreach(i=1:120)%dopar%{
      if(length(which(list.files(path = paste(l,"/",sep=""))==paste("All_Mapper_umap_p_",i,".RData",sep="")))==0){
        library(TDAmapper)
        library(igraph)
        library(umap)
        library(proxyC)
        library(cluster)
        library(bioDist)
        library(dbscan)
        lsmi_from_lsfi <- function( lsfi, num_intervals ) {
          # inputs:
          # lsfi = an integer in the range 1:prod(v)
          # num_intervals = c(i1,i1,...) a vector of numbers of intervals
          # output:
          # f+1 = a vector of multiindices with length filter_output_dim
          j <- c(1,num_intervals) # put 1 in front to make indexing easier in the product prod(j[1:k])
          f <- c()
          for (k in 1:length(num_intervals)) {
            # use lsfi-1 to shift from 1-based indexing to 0-based indexing
            f[k] <- floor( (lsfi-1) / prod(j[1:k])) %% num_intervals[k]
          }
          #print(f+1)
          # lsmi = f+1 = level set multi index
          return(f+1) # shift from 0-based indexing back to 1-based indexing
        }
        
        
        lsfi_from_lsmi <- function( lsmi, num_intervals ) {
          lsfi <- lsmi[1]
          if (length(num_intervals) > 1) {
            for (i in 2:length(num_intervals)) {
              lsfi <- lsfi + prod(num_intervals[1:(i-1)]) * (lsmi[i]-1)
            }
          }
          return(lsfi)
        }
        
        
        cluster_cutoff_at_first_empty_bin <- function(heights, diam, num_bins_when_clustering) {
          # if there are only two points (one height value), then we have a single cluster
          if (length(heights) == 1) {
            if (heights == diam) {
              cutoff <- Inf
              return(cutoff)
            }
          }
          
          bin_breaks <- seq(from=min(heights), to=diam,
                            by=(diam - min(heights))/num_bins_when_clustering)
          if (length(bin_breaks) == 1) { bin_breaks <- 1 }
          
          myhist <- hist(c(heights,diam), breaks=bin_breaks, plot=FALSE)
          z <- (myhist$counts == 0)
          if (sum(z) == 0) {
            cutoff <- Inf
            return(cutoff)
          } else {
            #  which returns the indices of the logical vector (z == TRUE), min gives the smallest index
            cutoff <- myhist$mids[ min(which(z == TRUE)) ]
            return(cutoff)
          }
          
        }
        cut_tree <- function(hcl, eps, core_dist){
          cuts <- unname(cutree(hcl, h=eps))
          #cuts[which(core_dist > eps)] <- 0 # Use core distance to distinguish noise
          return(cuts)
        }
        
        mapper=function (dist_object, filter_values, num_intervals, percent_overlap,
                         num_bins_when_clustering,clust_method="hclust")
        {
          filter_values <- data.frame(filter_values)
          num_points <- dim(filter_values)[1]
          filter_output_dim <- dim(filter_values)[2]
          num_levelsets <- prod(num_intervals)
          filter_min <- as.vector(sapply(filter_values, min))
          filter_max <- as.vector(sapply(filter_values, max))
          interval_width <- (filter_max - filter_min)/num_intervals
          vertex_index <- 0
          level_of_vertex <- c()
          points_in_vertex <- list()
          points_in_level_set <- vector("list", num_levelsets)
          vertices_in_level_set <- vector("list", num_levelsets)
          for (lsfi in 1:num_levelsets) {
            lsmi <- lsmi_from_lsfi(lsfi, num_intervals)
            lsfmin <- filter_min + (lsmi - 1) * interval_width -
              0.5 * interval_width * percent_overlap/100
            lsfmax <- lsfmin + interval_width + interval_width *
              percent_overlap/100
            for (point_index in 1:num_points) {
              if (all(lsfmin <= filter_values[point_index, ] &
                      filter_values[point_index, ] <= lsfmax)) {
                points_in_level_set[[lsfi]] <- c(points_in_level_set[[lsfi]],
                                                 point_index)
              }
            }
            points_in_this_level <- points_in_level_set[[lsfi]]
            num_points_in_this_level <- length(points_in_level_set[[lsfi]])
            if (num_points_in_this_level == 0) {
              num_vertices_in_this_level <- 0
            }
            if (num_points_in_this_level == 1) {
              num_vertices_in_this_level <- 1
              level_internal_indices <- c(1)
              level_external_indices <- points_in_level_set[[lsfi]]
            }
            if (num_points_in_this_level > 1) {
              if(clust_method=="hclust"){
                level_dist_object <- as.dist(as.matrix(dist_object)[points_in_this_level,
                                                                    points_in_this_level])
                level_max_dist <- max(level_dist_object)
                level_hclust <- hclust(level_dist_object, method = "single")
                level_heights <- level_hclust$height
                level_cutoff <- cluster_cutoff_at_first_empty_bin(level_heights,
                                                                  level_max_dist, num_bins_when_clustering)
                level_external_indices <- points_in_this_level[level_hclust$order]
                level_internal_indices <- as.vector(cutree(list(merge = level_hclust$merge,
                                                                height = level_hclust$height, labels = level_external_indices),
                                                           h = level_cutoff))
                level_external_indices <- points_in_this_level[sort(level_hclust$order)]
                num_vertices_in_this_level <- max(level_internal_indices)
              }
              if(clust_method=="hdbscan1"){
                level_dist_object <- as.dist(as.matrix(dist_object)[points_in_this_level,
                                                                    points_in_this_level])
                level_max_dist <- max(level_dist_object)
                level_hdbscan=hdbscan(level_dist_object,minPts = 2)
                level_heights <- level_hdbscan$hc$height
                level_cutoff <- cluster_cutoff_at_first_empty_bin(level_heights,
                                                                  level_max_dist, num_bins_when_clustering)
                level_external_indices <- points_in_this_level[level_hdbscan$hc$order]
                level_internal_indices <- as.vector(cut_tree(hcl = level_hdbscan$hc,eps = level_cutoff,core_dist = level_hdbscan$coredist))
                level_external_indices <- points_in_this_level[sort(level_hdbscan$hc$order)]
                num_vertices_in_this_level <- max(level_internal_indices)
                
              }
              if(clust_method=="hdbscan2"){
                level_dist_object <- as.dist(as.matrix(dist_object)[points_in_this_level,
                                                                    points_in_this_level])
                level_max_dist <- max(level_dist_object)
                level_hdbscan=hdbscan(level_dist_object,minPts = 2)
                level_heights <- level_hdbscan$hc$height
                level_cutoff <- cluster_cutoff_at_first_empty_bin(level_heights,
                                                                  level_max_dist, num_bins_when_clustering)
                level_external_indices <- points_in_this_level[sort(level_hdbscan$hc$order)]
                zeros=(which(level_hdbscan$cluster==0))
                zero_values=(max(level_hdbscan$cluster)+1):(max(level_hdbscan$cluster)+length(zeros))
                level_hdbscan$cluster[zeros]=zero_values
                level_internal_indices=level_hdbscan$cluster
                num_vertices_in_this_level <- max(level_internal_indices)
                
              }
              if(clust_method=="GN"){
                level_dist_object <- as.matrix(dist_object)[points_in_this_level,
                                                            points_in_this_level]
                g=graph.adjacency(level_dist_object,mode = "undirected",weighted = TRUE)
                level_GN=cluster_edge_betweenness(g,directed = FALSE)
                level_external_indices <- points_in_this_level
                level_internal_indices=level_GN$membership
                num_vertices_in_this_level <- max(level_internal_indices)
              }
              if(clust_method=="louvain"){
                level_dist_object <- as.matrix(dist_object)[points_in_this_level,
                                                            points_in_this_level]
                g=graph.adjacency(level_dist_object,mode = "undirected",weighted = TRUE)
                level_louvain=cluster_louvain(g,resolution = 0.75)
                level_external_indices <- points_in_this_level
                level_internal_indices=level_louvain$membership
                num_vertices_in_this_level <- max(level_internal_indices)
              }
              if(clust_method=="label"){
                level_dist_object <- as.matrix(dist_object)[points_in_this_level,
                                                            points_in_this_level]
                g=graph.adjacency(level_dist_object,mode = "undirected",weighted = TRUE)
                level_label=cluster_label_prop(g)
                level_external_indices <- points_in_this_level
                level_internal_indices=level_label$membership
                num_vertices_in_this_level <- max(level_internal_indices)
              }
              if(clust_method=="comp"){
                level_dist_object <- as.matrix(dist_object)[points_in_this_level,
                                                            points_in_this_level]
                g=graph.adjacency(level_dist_object,mode = "undirected",weighted = TRUE)
                level_comp=components(g)
                level_external_indices <- points_in_this_level
                level_internal_indices=level_comp$membership
                num_vertices_in_this_level <- max(level_internal_indices)
              }
            }
            if (num_vertices_in_this_level > 0) {
              vertices_in_level_set[[lsfi]] <- vertex_index + (1:num_vertices_in_this_level)
              for (j in 1:num_vertices_in_this_level) {
                vertex_index <- vertex_index + 1
                level_of_vertex[vertex_index] <- lsfi
                points_in_vertex[[vertex_index]] <- level_external_indices[level_internal_indices ==
                                                                             j]
              }
            }
          }
          
          adja <- mat.or.vec(vertex_index, vertex_index)
          for (lsfi in 1:num_levelsets) {
            lsmi <- lsmi_from_lsfi(lsfi, num_intervals)
            for (k in 1:filter_output_dim) {
              if (lsmi[k] < num_intervals[k]) {
                lsmi_adjacent <- lsmi + diag(filter_output_dim)[,
                                                                k]
                lsfi_adjacent <- lsfi_from_lsmi(lsmi_adjacent,
                                                num_intervals)
              }
              else {
                next
              }
              if (length(vertices_in_level_set[[lsfi]]) < 1 | length(vertices_in_level_set[[lsfi_adjacent]]) <
                  1) {
                next
              }
              for (v1 in vertices_in_level_set[[lsfi]]) {
                for (v2 in vertices_in_level_set[[lsfi_adjacent]]) {
                  adja[v1, v2] <- (length(intersect(points_in_vertex[[v1]],
                                                    points_in_vertex[[v2]])) > 0)
                  adja[v2, v1] <- adja[v1, v2]
                }
              }
            }
          }
          mapperoutput <- list(adjacency = adja, num_vertices = vertex_index,
                               level_of_vertex = level_of_vertex, points_in_vertex = points_in_vertex,
                               points_in_level_set = points_in_level_set, vertices_in_level_set = vertices_in_level_set)
          class(mapperoutput) <- "TDAmapper"
          return(mapperoutput)
        }
        mapper3D <- function(s,ppi_dist_matrix,emb_ppi,interv=10,per=20,bins=10,method="hclust"){
          ppi_mapper <- mapper(
            dist_object= ppi_dist_matrix,
            filter_values =list(emb_ppi[,1],emb_ppi[,2]),
            num_intervals =c(interv,interv),
            percent_overlap =per,
            num_bins_when_clustering = bins,
            clust_method = method)
          First.Example.mapper <- ppi_mapper
          
          MapperNodes <- mapperVertices(First.Example.mapper, (rownames(ppi_dist_matrix)) )
          return(MapperNodes)
        }
        print(paste("timestep: ",l," Interval: ",parameters[[i]][["interv"]]," Overlap: ",parameters[[i]][["per"]]," Current: ",i,"/120",sep=""))
        e=mapper3D(dat,ppi_dist_matrix,emb_ppi,parameters[[i]][["interv"]],parameters[[i]][["per"]],method="label")
        save(e,file=paste(l,"/All_Mapper_umap_p_",i,".RData",sep=""))
      }
    }
    stopCluster(clustery)
  }
}

Mapper()